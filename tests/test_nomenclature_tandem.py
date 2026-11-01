import random
import tempfile
from pathlib import Path
import unittest

import pysam
import test_structure as support
from asmsv_next.config import SVConfig
from asmsv_next.presets import resolve_preset
from asmsv_next.sequence_check import reverse_complement
from asmsv_next.structural_detection import detect_structures
from asmsv_next.genotype import cluster_and_genotype,cluster_structural_calls
from asmsv_next.structures import SVCandidate
from asmsv_next.nomenclature import duplication_metadata
from asmsv_next.vcf import write_vcf


class TandemNomenclatureTests(unittest.TestCase):
    def setUp(self):
        SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_SV_SIZE=30
        SVConfig.MIN_MAPQ=10;SVConfig.MERGE_INDEL_TOL=3;SVConfig.DETECTOR='structured'
        SVConfig.COPY_SOURCE_SEARCH='off';SVConfig.REPLACEMENT_REALIGN='off'
        rng=random.Random(885);self.ref=''.join(rng.choice('ACGT') for _ in range(15000))

    def tandem(self,length,filler='',reverse=False):
        start=2000;end=start+length
        q=self.ref[:end]+filler+self.ref[start:end+2000]
        chain=[support.block(0,[(7,end)],0,'+',len(q)),
               support.block(end+len(filler),[(7,length+2000)],start,'+',len(q))]
        if reverse:
            chain=[support.block(len(q)-b.qe,[(op,n) for op,n in b.cigar if op not in (4,5)],
                b.ref_start,'-',len(q)) for b in reversed(chain)]
            q=reverse_complement(q)
        return q,detect_structures(chain,'-' if reverse else '+',{0,1},'q',q,support.Reference(self.ref))

    def test_tandem_type_has_no_1000bp_switch_and_projection_spells_allele(self):
        for length in (30,99,500,999,1000,1200):
            for reverse in (False,True):
                q,result=self.tandem(length,reverse=reverse)
                self.assertEqual([c.type for c in result.calls],['DUP'])
                parent=result.calls[0]
                self.assertEqual((parent.info['DUP_MODE'],parent.info['DUP_ORIENT']),('TANDEM','DIRECT'))
                self.assertEqual(parent.sv_len,length);self.assertEqual(result.cigar_scans,2)
                projection=result.projections[0];target=projection.pos
                reconstructed=self.ref[:target]+projection.info['SEQ']+self.ref[target:4000+length]
                self.assertEqual(reconstructed,reverse_complement(q) if reverse else q)

    def test_nonpositive_legacy_gain_does_not_create_a_duplication(self):
        q=self.ref[:5000]
        chain=[support.block(0,[(7,2500)],0,'+',len(q)),support.block(2420,[(7,2000)],2460,'+',len(q))]
        result=detect_structures(chain,'+',{0,1},'q',q,support.Reference(self.ref))
        self.assertFalse(result.calls);self.assertFalse(result.projections)
        self.assertIn('QUERY_OWNERSHIP_CONFLICT',result.contexts[0]['issues'])
        self.assertEqual(result.cigar_scans,2)

    def test_positive_gap_can_be_an_exact_anchor_extension(self):
        for filler,method in [(self.ref[1992:2000],'RIGHT_ANCHOR_EXTENSION'),
                              (self.ref[2500:2508],'LEFT_ANCHOR_EXTENSION')]:
            q,result=self.tandem(500,filler)
            self.assertEqual([c.type for c in result.calls],['DUP'])
            self.assertEqual(result.calls[0].info['DUP_NORMALIZATION'],method)
            projection=result.projections[0];target=projection.pos
            self.assertEqual(self.ref[:target]+projection.info['SEQ']+self.ref[target:4500],q)
            self.assertTrue(result.contexts[0]['boundary_extensions'])

    def test_novel_junction_is_explicit_and_does_not_claim_pure_tandem(self):
        for length in (20,30,50):
            for reverse in (False,True):
                q,result=self.tandem(500,'G'*length,reverse)
                parent=next(c for c in result.calls if c.type in ('DUP','CPX'))
                self.assertEqual(parent.info['DUP_MODE'],'ADJACENT')
                self.assertEqual(parent.info['JUNCTION_INS_LEN'],length)
                if length>=30:
                    self.assertEqual(parent.type,'CPX');self.assertEqual(parent.info['CPX_TYPE'],'COMPLEX')
                    self.assertEqual(parent.info['CPX_COMPONENTS'],'DUP,INS')
                else:self.assertEqual(parent.type,'DUP')
                projection=result.projections[0];target=projection.pos
                self.assertEqual(self.ref[:target]+projection.info['SEQ']+self.ref[target:4500],reverse_complement(q) if reverse else q)

    def test_edit_in_left_or_right_duplicate_has_distinct_target_and_allele(self):
        left_query=self.ref[:2500]+self.ref[2700:3500]+self.ref[2000:5500]
        right_query=self.ref[:3500]+self.ref[2000:2500]+self.ref[2700:5500]
        left_chain=[support.block(0,[(7,2500),(2,200),(7,800)],0,'+',len(left_query)),
                    support.block(3300,[(7,3500)],2000,'+',len(left_query))]
        right_chain=[support.block(0,[(7,3500)],0,'+',len(right_query)),
                     support.block(3500,[(7,500),(2,200),(7,2800)],2000,'+',len(right_query))]
        results=[detect_structures(chain,'+',{0,1},'q',q,support.Reference(self.ref))
                 for chain,q in [(left_chain,left_query),(right_chain,right_query)]]
        parents=[next(c for c in r.calls if c.type=='CPX') for r in results]
        self.assertEqual([p.info['TARGET_POS'] for p in parents],[2000,3500])
        self.assertNotEqual(parents[0].info['STRUCTURE_KEY'],parents[1].info['STRUCTURE_KEY'])
        self.assertEqual(len(cluster_and_genotype([parents[0]],[parents[1]])),2)
        for result,q in zip(results,[left_query,right_query]):
            projection=result.projections[0];p=projection.pos
            self.assertEqual(self.ref[:p]+projection.info['SEQ']+self.ref[p:5500],q)

    def test_cpx_names_follow_reference_flanks_not_query_orientation(self):
        q=self.ref[:1900]+reverse_complement(self.ref[2000:3000])+self.ref[3000:5000]
        chain=[support.block(0,[(7,1900)],0,'+',len(q)),support.block(1900,[(7,1000)],2000,'-',len(q)),
               support.block(2900,[(7,2000)],3000,'+',len(q))]
        first=detect_structures(chain,'+',{0,2},'q',q,support.Reference(self.ref))
        reversed_chain=[support.block(len(q)-b.qe,[(op,n) for op,n in b.cigar if op not in (4,5)],
            b.ref_start,'+' if b.strand=='-' else '-',len(q)) for b in reversed(chain)]
        second=detect_structures(reversed_chain,'-',{0,2},'q',reverse_complement(q),support.Reference(self.ref))
        parents=[next(c for c in r.calls if c.type=='CPX') for r in (first,second)]
        self.assertTrue(all(p.info['CPX_TYPE']=='delINV' for p in parents))
        self.assertEqual(parents[0].info['STRUCTURE_KEY'],parents[1].info['STRUCTURE_KEY'])
        merged=cluster_and_genotype([parents[0]],[parents[1]])
        self.assertEqual(len(merged),1);self.assertEqual(merged[0].genotype,'1|1')

    def test_adjacent_inverted_duplication_is_not_direct_tandem(self):
        for target,subtype in [(3000,'piDUP_FR'),(2000,'piDUP_RF')]:
            info=duplication_metadata('chr1',2000,3000,'chr1',target,'INVERTED')
            self.assertEqual(info['DUP_MODE'],'ADJACENT');self.assertEqual(info['SUBTYPE'],subtype)

    def test_vcf_exposes_component_terms_without_old_type_labels(self):
        q,result=self.tandem(500,'G'*50)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'calls.vcf'
            write_vcf(result.calls,str(path),{'chr1':len(self.ref)},support.Reference(self.ref))
            text=path.read_text()
            self.assertNotIn('CSV_TYPE',text);self.assertNotIn('COPY_',text);self.assertNotIn('TEMPLATED_CHAIN',text)
            rows=list(pysam.VariantFile(str(path)));parent=next(r for r in rows if r.info['SVTYPE']=='CPX')
            self.assertEqual(parent.info['CPX_COMPONENTS'],('DUP','INS'))

    def test_added_duplicate_deletion_does_not_genotype_as_reference_deletion(self):
        parent=SVCandidate('CPX','chr1',2000,900,'h1',dict(ROLE='EVENT',EVENT='p1',CPX_TYPE='COMPLEX',CPX_COMPONENTS='DUP,DEL'))
        parent.id='p1'
        child=SVCandidate('DEL','chr1',2400,100,'h1',dict(ROLE='COMPONENT',PARENT_EVENT='p1',ALLELE_CONTEXT='ADDED_DUPLICATE'))
        simple=SVCandidate('DEL','chr1',2400,100,'h2',dict(ROLE='SIMPLE'))
        results,aliases=cluster_structural_calls([parent,child],[simple])
        dels=[c for c in results if c.type=='DEL']
        self.assertEqual(len(dels),2)
        self.assertEqual({(c.info['ROLE'],c.genotype) for c in dels},{('COMPONENT','1|0'),('SIMPLE','0|1')})

    def test_cigar_and_split_tandem_share_native_type_without_changing_ins_atom(self):
        SVConfig.COPY_SOURCE_SEARCH='adjacent'
        q,split=self.tandem(500)
        chain=[support.block(0,[(7,2500),(1,500),(7,2000)],0,'+',len(q))]
        cigar=detect_structures(chain,'+',{0},'q',q,support.Reference(self.ref))
        native=[c for c in cigar.calls if c.info['ROLE']!='COMPONENT']
        self.assertEqual([c.type for c in native],['DUP'])
        raw=next(c for c in cigar.calls if c.type=='INS')
        self.assertEqual(raw.info['SEQ'],self.ref[2000:2500]);self.assertEqual(raw.sv_len,500)
        self.assertEqual(raw.info['ALLELE_CONTEXT'],'REFERENCE')
        self.assertEqual(cigar.cigar_scans,1);self.assertFalse(cigar.projections)
        merged=cluster_and_genotype([native[0]],[split.calls[0]])
        self.assertEqual(len(merged),1);self.assertEqual(merged[0].genotype,'1|1')
        gap_chain=[support.block(0,[(7,2500)],0,'+',len(q)),support.block(3000,[(7,2000)],2500,'+',len(q))]
        gap=detect_structures(gap_chain,'+',{0,1},'q',q,support.Reference(self.ref))
        gap_native=[c for c in gap.calls if c.info['ROLE']!='COMPONENT']
        self.assertEqual([c.type for c in gap_native],['DUP'])
        self.assertFalse(any('OPEN' in reason for ctx in gap.contexts for reason in ctx['issues']))
        self.assertEqual(len(cluster_and_genotype([native[0]],gap_native)),1)

    def test_multiple_cigar_duplications_in_one_block_have_distinct_contexts(self):
        SVConfig.COPY_SOURCE_SEARCH='adjacent'
        q=self.ref[:2500]+self.ref[2000:2500]+self.ref[2500:4000]+self.ref[3500:4000]+self.ref[4000:6000]
        chain=[support.block(0,[(7,2500),(1,500),(7,1500),(1,500),(7,2000)],0,'+',len(q))]
        result=detect_structures(chain,'+',{0},'q',q,support.Reference(self.ref))
        self.assertEqual(len({c['context_id'] for c in result.contexts}),len(result.contexts))
        parents=[c for c in result.calls if c.type=='DUP']
        self.assertEqual(len(parents),2);self.assertEqual(len({p.id for p in parents}),2)


if __name__=='__main__':unittest.main()
