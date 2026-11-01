import tempfile
from pathlib import Path
import unittest
import pysam

from asmsv_next.config import SVConfig
from asmsv_next.structures import SVCandidate
from asmsv_next.vcf import write_vcf
from asmsv_next.vcf_annotations import compact_record


class AnnotationTests(unittest.TestCase):
    def setUp(self):SVConfig.DETECTOR='structured'

    def render(self,calls,view='main',level='minimal',sink=None):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'calls.vcf'
            write_vcf(calls,str(path),{'chr1':10000},None,view=view,info_level=level,record_sink=sink)
            text=path.read_text();rows=list(pysam.VariantFile(str(path)))
            return text,rows

    def test_tandem_does_not_repeat_alt_type_and_algorithm_state(self):
        c=SVCandidate('DUP','chr1',100,500,'q',dict(ROLE='EVENT',EVENT='p',SUBTYPE='DUP:TANDEM',DUP_MODE='TANDEM',
            DUP_ORIENT='DIRECT',SOURCE='STRUCTURE',SOURCE_STATE='RETAINED',FRAME='DUP',NET_LEN=500,
            DUP_TEMPLATE_LEN=500,ATOM_ID='a',CONTEXT_ID='ctx',RESOLUTION='RESOLVED'))
        c.id='p'
        text,rows=self.render([c]);info=dict(rows[0].info)
        self.assertEqual(rows[0].alts,('<DUP:TANDEM>',))
        self.assertNotIn('DUP_TYPE',info);self.assertNotIn('DUP_MODE',info);self.assertNotIn('NET_LEN',info)
        self.assertNotIn('ATOM_ID',info);self.assertNotIn('EVENT',info);self.assertNotIn('RESOLUTION',info)

    def test_copy_local_component_keeps_context_and_standard_event_link(self):
        c=SVCandidate('DEL','chr1',200,100,'q',dict(ROLE='COMPONENT',PARENT_EVENT='p',
            ALLELE_CONTEXT='ADDED_DUPLICATE',ATOM_ID='a',CONTEXT_ID='ctx'))
        captured=[];text,rows=self.render([c],view='components',sink=captured.append)
        self.assertEqual(rows[0].info['EVENT'],'p');self.assertEqual(rows[0].info['EDIT_CONTEXT'],'ADDED_DUPLICATE')
        self.assertNotIn('PARENT_EVENT',dict(rows[0].info));self.assertNotIn('ALLELE_CONTEXT',text)
        self.assertEqual(captured[0]['record_id'],rows[0].id);self.assertEqual(captured[0]['info']['ATOM_ID'],'a')

    def test_distinct_haplotype_parents_are_not_forced_into_one_event(self):
        c=SVCandidate('INS','chr1',200,100,'q',dict(ROLE='PROJECTION',PARENT_EVENT='p1,p2',SEQ='A'*100))
        text,rows=self.render([c],view='indels');info=dict(rows[0].info)
        self.assertEqual(info['PARENT_EVENT'],('p1','p2'));self.assertNotIn('EVENT',info)

    def test_generic_complex_label_is_omitted_but_component_edits_remain(self):
        c=SVCandidate('CPX','chr1',100,900,'q',dict(ROLE='EVENT',CPX_TYPE='COMPLEX',CPX_COMPONENTS='DUP,DEL',
            CPX_INTERVALS='DUP_chr1:101-1100,DEL_chr1:201-300',INS_LEN=0,DEL_LEN=100,NET_LEN=900),end=1100)
        text,rows=self.render([c]);info=dict(rows[0].info)
        self.assertNotIn('CPX_TYPE',info);self.assertEqual(info['CPX_COMPONENTS'],('DUP','DEL'))
        self.assertEqual(info['EDIT_DEL_LEN'],100);self.assertNotIn('EDIT_INS_LEN',info)

    def test_full_profile_has_diagnostics_and_all_core_columns_are_identical(self):
        c=SVCandidate('INS','chr1',200,100,'q',dict(SOURCE='CIGAR',ROLE='SIMPLE',ATOM_ID='a',SEQ='A'*100))
        _,minimal=self.render([c]);_,full=self.render([c],level='full')
        self.assertEqual(full[0].info['DETECTION_METHOD'],'CIGAR');self.assertEqual(full[0].info['ATOM_ID'],('a',))
        for row in (minimal[0],full[0]):
            self.assertEqual((row.chrom,row.pos,row.ref,row.alts,row.samples[0]['GT']),('chr1',200,'N',('N'+'A'*100,),(None,None)))

    def test_breakend_mates_and_alt_are_kept_without_duplicate_ports(self):
        c=SVCandidate('BND','chr1',200,0,'q',dict(CHR2='chr1',POS2=900,STRANDS='++',
            ROLE='JUNCTION',JUNCTION_ID='j',RESOLUTION='PARTIAL'),filter_str='PARTIAL')
        c.id='j.a';c.mateid='j.b'
        _,rows=self.render([c],view='junctions');info=dict(rows[0].info)
        self.assertEqual(info['MATEID'],('j.b',));self.assertIn('chr1:900',rows[0].alts[0])
        self.assertNotIn('CHR2',info);self.assertNotIn('STRANDS',info);self.assertNotIn('RESOLUTION',info)

    def test_known_source_orientation_and_unknown_retention_are_not_discarded(self):
        c=SVCandidate('INS','chr1',900,100,'q',dict(INS_ORIGIN='TEMPLATED',SOURCE_STATE='UNKNOWN',
            SOURCE_CHR='chr1',SOURCE_START=100,SOURCE_END=200,SOURCE_ORIENT='INVERTED',SEQ='A'*100))
        _,rows=self.render([c]);info=dict(rows[0].info)
        self.assertEqual(info['SOURCE_ORIENT'],'INVERTED');self.assertEqual(info['SOURCE_STATE'],'UNKNOWN')

    def test_copy_origin_of_an_inversion_insertion_component_is_preserved(self):
        c=SVCandidate('INS','chr1',900,100,'q',dict(ROLE='COMPONENT',PARENT_EVENT='inv',
            ALLELE_CONTEXT='REFERENCE',INS_ORIGIN='DUPLICATION',SOURCE_STATE='RETAINED',
            SOURCE_CHR='chr1',SOURCE_START=100,SOURCE_END=200,DUP_ORIENT='INVERTED',SEQ='A'*100))
        _,rows=self.render([c],view='components');info=dict(rows[0].info)
        self.assertEqual(info['SOURCE_START'],100);self.assertEqual(info['SOURCE_END'],200)
        self.assertEqual(info['SOURCE_ORIENT'],'INVERTED');self.assertEqual(info['INS_ORIGIN'],'DUPLICATION')
