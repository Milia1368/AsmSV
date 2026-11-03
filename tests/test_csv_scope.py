import json
from pathlib import Path
import tempfile
import unittest

import pysam

from asmsv_next.csv_scope import csv_class
from asmsv_next.cpx_export import load_events, main, view_payload
from asmsv_next.cpx_diagram import diagrams_for_event


HEADER='''##fileformat=VCFv4.2
##contig=<ID=chr2>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=DUP_TYPE,Number=1,Type=String,Description="Subtype">
##INFO=<ID=DUP_ORIENT,Number=1,Type=String,Description="Direction">
##INFO=<ID=TARGET_CHR,Number=1,Type=String,Description="Target chr">
##INFO=<ID=TARGET_POS,Number=1,Type=Integer,Description="Target cut">
##INFO=<ID=INS_LEN,Number=1,Type=Integer,Description="Inserted length">
##INFO=<ID=DEL_LEN,Number=1,Type=Integer,Description="Deleted length">
##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Intervals">
##INFO=<ID=CPX_COMPONENTS,Number=.,Type=String,Description="Components">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
'''
USER_ROW='chr2\t139197368\tEVENT.32f0f792b4099c540029b407\tG\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=139347370;SVLEN=150002;TARGET_CHR=chr2;TARGET_POS=53343523;DUP_ORIENT=INVERTED;DUP_TYPE=dDUP\tGT:PS\t1|0:10001\n'


class CsvScopeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.vcf=self.root/'calls.vcf'

    def test_current_production_classes_are_explicit(self):
        for info,expected in [({'SVTYPE':'CPX'},'COMPOSITE'),({'SVTYPE':'DELINS'},'REPLACEMENT'),
          ({'SVTYPE':'DUP','DUP_TYPE':'piDUP_FR'},'INVERTED_DUPLICATION'),
          ({'SVTYPE':'DUP','DUP_TYPE':'piDUP_RF'},'INVERTED_DUPLICATION'),
          ({'SVTYPE':'DUP','DUP_TYPE':'dDUP','DUP_ORIENT':'INVERTED'},'INVERTED_DISPERSED_DUPLICATION'),
          ({'SVTYPE':'DUP','DUP_TYPE':'dDUP','DUP_ORIENT':'DIRECT'},'DISPERSED_DUPLICATION'),
          ({'SVTYPE':'DUP','DUP_ORIENT':'INVERTED'},'INVERTED_DUPLICATION')]:
            self.assertEqual(csv_class(info),expected)

    def test_cpx_class_describes_structure_instead_of_repeating_svtype(self):
        self.assertEqual(csv_class(dict(SVTYPE='CPX',CPX_COMPONENTS='INV,DEL,DEL')),'INV_WITH_EDITS')
        self.assertEqual(csv_class(dict(SVTYPE='CPX',CPX_COMPONENTS='DUP,DEL,INS',
                                       CPX_INTERVALS='DUP_chr2:101-500,DEL_chr2:201-250')),'DUP_WITH_EDITS')
        self.assertEqual(csv_class(dict(SVTYPE='CPX',CPX_COMPONENTS='DUP,INS',INS_ORIGIN='TEMPLATED',
                                       CPX_INTERVALS='DUP_chr2:101-500,INS_chr3:201-250')),'MULTISOURCE_INSERTION')
        # The same type list alone does not establish a two-source layout.
        self.assertEqual(csv_class(dict(SVTYPE='CPX',CPX_COMPONENTS='DUP,INS')),'COMPOSITE')

    def test_basic_and_unresolved_source_annotations_are_not_promoted(self):
        for kind in ['INS','DEL','INV','DUP','TRA','BND','UNKNOWN']:
            self.assertIsNone(csv_class(dict(SVTYPE=kind)))
        self.assertIsNone(csv_class(dict(SVTYPE='INS',INS_ORIGIN='TEMPLATED',SOURCE_ORIENT='INVERTED')))
        self.assertIsNone(csv_class(dict(SVTYPE='DUP'),'<DUP:TANDEM>'))
        self.assertIsNone(csv_class(dict(SVTYPE='DUP',DUP_TYPE='DUP:TANDEM')))

    def test_user_native_inverted_disperse_copy_is_exported_without_relabeling(self):
        self.vcf.write_text(HEADER+USER_ROW)
        prefix=self.root/'calls.csv'
        self.assertEqual(main(['--vcf',str(self.vcf),'--out-prefix',str(prefix)]),0)
        with pysam.VariantFile(str(prefix)+'.verix.vcf') as stream:
            rows=list(stream)
        self.assertEqual(len(rows),1);r=rows[0]
        self.assertEqual(r.info['SVTYPE'],'DUP')
        self.assertEqual(r.info['ASMSV_CSV_CLASS'],'INVERTED_DISPERSED_DUPLICATION')
        self.assertEqual(r.info['BKPS'],('DUP:150002-chr2:139197368-chr2:139347370','TARGET:0-chr2:53343523'))
        self.assertEqual(r.samples['S']['GT'],(1,0));self.assertTrue(r.samples['S'].phased)
        self.assertEqual(r.samples['S']['PS'],10001)
        items=[json.loads(l) for l in Path(str(prefix)+'.jsonl').read_text().splitlines()]
        self.assertEqual(items[0]['scope'],'csv')
        self.assertEqual(items[1]['record_kind'],'CSV_EVENT')
        m=diagrams_for_event(items[1])[0]
        self.assertEqual((m['status'],m['layout'],m['copy_direction']),('supported','source_target','-'))

    def test_strict_scope_excludes_native_copy(self):
        self.vcf.write_text(HEADER+USER_ROW)
        self.assertEqual(load_events(self.vcf,None,scope='cpx'),[])

    def test_replacement_uses_both_lengths_even_when_net_svlen_is_one(self):
        row='chr2\t100\treplacement\tN\t<DELINS>\t.\tPASS\tSVTYPE=DELINS;END=210;DEL_LEN=110;INS_LEN=111;SVLEN=1\tGT\t1/.\n'
        self.vcf.write_text(HEADER+row)
        events=load_events(self.vcf,None)
        self.assertEqual(events[0]['bkps'],['DEL:110-chr2:100-chr2:210','INS:111-chr2:100'])
        model=diagrams_for_event(events[0])[0]
        self.assertEqual(model['status'],'supported')
        self.assertEqual([s['length'] for s in model['summaries']],[110,111])
        self.assertEqual(events[0]['main']['info']['SVLEN'],'1')

    def test_pi_dup_direction_is_known_even_when_minimal_profile_omits_it(self):
        row='chr2\t100\tpi\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=200;SVLEN=100;DUP_TYPE=piDUP_RF;TARGET_CHR=chr2;TARGET_POS=100\tGT\t0/1\n'
        self.vcf.write_text(HEADER+row)
        events=load_events(self.vcf,None)
        model=diagrams_for_event(events[0])[0]
        self.assertEqual(model['status'],'supported')
        self.assertEqual([g['role'] for g in model['alt']],['added','retained'])
        self.assertEqual(model['alt'][0]['direction'],'-')
        self.assertNotIn('DUP_ORIENT',events[0]['main']['info'])

    def test_copy_ins_basis_is_not_counted_as_novel_insertion(self):
        self.vcf.write_text(HEADER+USER_ROW)
        comps=self.root/'components.vcf'
        comps.write_text(HEADER+'chr2\t53343523\tbasis\tN\t<INS>\t.\tPASS\tSVTYPE=INS;END=53343523;SVLEN=150002;INS_ORIGIN=DUPLICATION;EVENT=EVENT.32f0f792b4099c540029b407\tGT\t1/.\n')
        event=load_events(self.vcf,comps)[0]
        self.assertEqual(len(event['components']),1)
        self.assertEqual(len(event['bkps']),2)
        self.assertFalse(any(e.startswith('INS:') for e in event['bkps']))

    def test_native_copy_evidence_joins_alias_and_preserves_projection(self):
        self.vcf.write_text(HEADER+USER_ROW)
        path=self.root/'evidence.jsonl'
        e=dict(event_id='internal',svtype='DUP',context_id='ctx',info={},components=[],haplotype='hap1')
        c=dict(context_id='ctx',insertion_projections=[dict(length=150002)])
        path.write_text(json.dumps(dict(record_kind='CONTIG',events=[e],contexts=[c],atoms=[]))+'\n'+
                        json.dumps(dict(record_kind='EVENT_ALIASES',aliases={'internal':'EVENT.32f0f792b4099c540029b407'})))
        events=load_events(self.vcf,None,path)
        self.assertEqual(len(events[0]['structures']),1)
        self.assertEqual(diagrams_for_event(events[0])[0]['projections'][0]['length'],150002)


if __name__=='__main__':unittest.main()
