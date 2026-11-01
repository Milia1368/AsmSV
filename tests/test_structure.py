import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import random
from types import SimpleNamespace
import unittest

from asmsv_next.alignment import AlignedSegment
from asmsv_next.config import SVConfig
from asmsv_next.presets import resolve_preset
from asmsv_next.structures import SVCandidate
from asmsv_next.ssv import make_gap_caller
from fixtures.legacy_detection import resolve_chain_svs
from asmsv_next.structural_detection import detect_structures
from asmsv_next.sequence_check import check_sequence,reverse_complement
from asmsv_next.output_views import junction_candidates
from asmsv_next.genotype import cluster_and_genotype

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_cigar',HERE/'fixtures/legacy_cigar.py')
oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)


class Reference:
    def __init__(self,sequence):self.sequence=sequence
    def fetch(self,chrom,start,end):return self.sequence[start:end]
    def get_reference_length(self,chrom):return len(self.sequence)


def block(qs,body,rs,strand,qlen):
    size=sum(n for op,n in body if op in (0,1,7,8));qe=qs+size
    left,right=(qs,qlen-qe) if strand=='+' else (qlen-qe,qs)
    cigar=([(4,left)] if left else [])+body+([(4,right)] if right else [])
    return AlignedSegment(dict(q_name='q',ref_name='chr1',ref_start=rs,
        ref_end=rs+sum(n for op,n in body if op in (0,2,7,8)),cigar=cigar,
        is_reverse=strand=='-',mapq=60,query_length=qlen))


class StructureTests(unittest.TestCase):
    def setUp(self):
        SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_SV_SIZE=30
        SVConfig.MIN_MAPQ=10;SVConfig.MERGE_INDEL_TOL=3
        SVConfig.COPY_SOURCE_SEARCH='off';SVConfig.REPLACEMENT_REALIGN='off'
        rng=random.Random(998);self.ref=''.join(rng.choice('ACGT') for _ in range(12000))

    def test_lifted_gap_kernel_is_exact_frozen_ast(self):
        expected=json.loads((HERE/'fixtures/gap_kernel_ast.json').read_text())
        source=(HERE.parent/'src/asmsv_next/ssv.py').read_text();tree=ast.parse(source)
        wrapper=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='make_gap_caller')
        for node in wrapper.body:
            if isinstance(node,ast.FunctionDef):
                # Python 3.12 adds empty type-parameter metadata. Canonicalize
                # only that empty field; retain actual generic declarations.
                for child in ast.walk(node):
                    if getattr(child,'type_params',None)==[]:
                        del child.type_params
                self.assertEqual(hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest(),expected[node.name])

    def test_gap_kernel_boundaries_and_reverse_strand_match_original(self):
        query='ACGT'*2500;kernel=make_gap_caller('q',query)
        signature=lambda c:None if c is None else (c.type,c.chrom,c.pos,c.sv_len,c.end,c.info,c.q_start,c.q_end,c.reference_gap)
        for strand in ('+','-'):
            for gq in (-101,-99,-2,0,29,30,99,100,400,1000):
                for gr in (-2000,-1000,-31,-30,-1,0,30,39,40,99,100,500,150000,150001):
                    a=SimpleNamespace(qs=0,qe=2000,ref_name='chr1',ref_start=500000,ref_end=502000,strand=strand,is_backbone=True)
                    rs,re=(502000+gr,504000+gr) if strand=='+' else (498000-gr,500000-gr)
                    b=SimpleNamespace(qs=2000+gq,qe=4000+gq,ref_name='chr1',ref_start=rs,ref_end=re,strand=strand,is_backbone=True)
                    old=resolve_chain_svs([a,b],strand,{0,1},'q',query)
                    new=kernel(a,b,'chr1',strand,'backbone')
                    self.assertEqual(signature(new),signature(old[0] if old else None))

    def test_cigar_oracle_and_provenance_randomized(self):
        rng=random.Random(44)
        for _ in range(200):
            body=[(7,100)]
            for j in range(12):body.extend([(rng.choice((1,2,8)),rng.choice((1,3,20,29,30,31,50,100))),(7,rng.choice((1,3,4,17)))])
            body.append((7,100))
            for reverse in (False,True):
                cigar=[(5,71),(4,13)]+body+[(4,19),(5,53)]
                data=dict(q_name='q',ref_name='chr1',ref_start=1000,cigar=cigar,mapq=60,is_reverse=reverse)
                data['ref_end']=1000+sum(n for op,n in cigar if op in (0,2,7,8))
                old=oracle.LegacySegment(data);new=AlignedSegment(data)
                sequence='ACGT'*1000
                expected=old.parse_intra_alignment_sv(sequence,30,3)
                self.assertEqual(new.parse_intra_alignment_sv(sequence,30,3),expected)
                annotated=new.parse_intra_alignment_sv(sequence,30,3,with_provenance=True)
                self.assertEqual([{k:v for k,v in x.items() if k!='evidence'} for x in annotated],expected)
                for atom in annotated:
                    for op,qs,qe,rs,re in atom['evidence']['parts']:
                        self.assertTrue(new.qs<=qs<=qe<=new.qe)

    def test_single_branch_inv_is_simple_and_scanned_once(self):
        q=self.ref[:2000]+reverse_complement(self.ref[2000:3000])+self.ref[3000:5000]
        chain=[block(0,[(7,2000)],0,'+',len(q)),block(2000,[(7,1000)],2000,'-',len(q)),block(3000,[(7,2000)],3000,'+',len(q))]
        r=detect_structures(chain,'+',{0,2},'q',q,Reference(self.ref))
        self.assertEqual([c.type for c in r.calls],['INV']);self.assertEqual(r.cigar_scans,3)
        self.assertEqual(r.calls[0].sv_len,1000);self.assertEqual(r.calls[0].info['NET_LEN'],0)

    def test_unaligned_inversion_source_is_not_complete_inv_delins(self):
        q=self.ref[:2000]+reverse_complement(self.ref[2000:3500])+self.ref[3500:5500]
        chain=[block(0,[(7,2000)],0,'+',len(q)),block(2800,[(7,700)],2000,'-',len(q)),
            block(3500,[(7,2000)],3500,'+',len(q))]
        result=detect_structures(chain,'+',{0,2},'q',q,Reference(self.ref))
        self.assertFalse(any(c.type=='CPX' for c in result.calls))
        self.assertEqual(sorted(c.type for c in result.calls),['DEL','INS','INV'])
        self.assertIn('INVERSION_BOUNDARY_SOURCE_UNRESOLVED',result.contexts[0]['issues'])
        self.assertEqual(result.cigar_scans,3)

    def test_no_backbone_retains_cigar_evidence_without_inventing_junction(self):
        q=self.ref[:2000]+self.ref[2100:4000]
        chain=[block(0,[(7,2000),(2,100),(7,1900)],0,'+',len(q))]
        result=detect_structures(chain,'+',set(),'q',q,Reference(self.ref))
        self.assertEqual([c.type for c in result.calls],['DEL'])
        self.assertEqual(result.cigar_scans,1)
        self.assertEqual(result.contexts[0]['issues'],['OPEN_LEFT','OPEN_RIGHT'])
        self.assertTrue(result.atoms[0]['contexts'])

    def test_inv_internal_deletion_has_one_parent_and_one_component(self):
        q=self.ref[:2000]+reverse_complement(self.ref[2700:3500])+reverse_complement(self.ref[2000:2500])+self.ref[3500:5500]
        chain=[block(0,[(7,2000)],0,'+',len(q)),block(2000,[(7,500),(2,200),(7,800)],2000,'-',len(q)),block(3300,[(7,2000)],3500,'+',len(q))]
        r=detect_structures(chain,'+',{0,2},'q',q,Reference(self.ref))
        self.assertEqual(len([c for c in r.calls if c.type=='CPX']),1)
        d=next(c for c in r.calls if c.type=='DEL');p=next(c for c in r.calls if c.type=='CPX')
        self.assertEqual((p.info['CPX_COMPONENTS'],d.sv_len,d.info['PARENT_EVENT']),('INV,DEL',200,p.id))
        self.assertEqual(r.cigar_scans,len(chain))

    def test_backbone_copy_edit_does_not_need_branch(self):
        q=self.ref[:3500]+self.ref[2000:2500]+self.ref[2700:3500]+self.ref[3500:5500]
        chain=[block(0,[(7,3500)],0,'+',len(q)),block(3500,[(7,500),(2,200),(7,2800)],2000,'+',len(q))]
        r=detect_structures(chain,'+',{0,1},'q',q,Reference(self.ref))
        parent=next(c for c in r.calls if c.type=='CPX')
        self.assertEqual(parent.info['CPX_COMPONENTS'],'DUP,DEL');self.assertEqual(parent.info['NET_LEN'],1300)
        self.assertEqual(r.cigar_scans,2)
        self.assertTrue(all(c.info.get('ROLE')=='COMPONENT' for c in r.calls if c.type in ('DUP','DEL')))

    def test_copy_retention_clips_query_overlap_at_long_run_boundary(self):
        q=self.ref[:3500]+self.ref[2000:2500]+self.ref[2700:3500]+self.ref[3500:5500]
        # Two exact boundary bases can be represented in both alignments.
        chain=[block(0,[(7,3502)],0,'+',len(q)),block(3500,[(7,500),(2,200),(7,2800)],2000,'+',len(q))]
        result=detect_structures(chain,'+',{0,1},'q',q,Reference(self.ref))
        parent=next(c for c in result.calls if c.type=='CPX')
        self.assertEqual(parent.info['CPX_COMPONENTS'],'DUP,DEL')
        self.assertEqual(result.cigar_scans,2)

    def test_retained_inverted_copy_is_simple_dup(self):
        q=self.ref[:3500]+reverse_complement(self.ref[2000:3000])+self.ref[3500:5500]
        chain=[block(0,[(7,3500)],0,'+',len(q)),block(3500,[(7,1000)],2000,'-',len(q)),block(4500,[(7,2000)],3500,'+',len(q))]
        r=detect_structures(chain,'+',{0,2},'q',q,Reference(self.ref))
        self.assertEqual([c.type for c in r.calls],['DUP']);self.assertEqual(r.calls[0].info['DUP_MODE'],'DISPERSED')

    def test_source_insertion_with_sv_sized_target_loss_is_not_simple_dup(self):
        q=self.ref[:3500]+reverse_complement(self.ref[2000:3000])+self.ref[3570:5570]
        chain=[block(0,[(7,3500)],0,'+',len(q)),block(3500,[(7,1000)],2000,'-',len(q)),
            block(4500,[(7,2000)],3570,'+',len(q))]
        result=detect_structures(chain,'+',{0,2},'q',q,Reference(self.ref))
        self.assertFalse(any(c.type in ('DUP','CPX') for c in result.calls))
        self.assertIn('TARGET_REPLACEMENT_WITH_SOURCE_PATH',result.contexts[0]['issues'])

    def test_whole_query_reverse_complement_preserves_copy_edit(self):
        q=self.ref[:3500]+self.ref[2000:2500]+self.ref[2700:3500]+self.ref[3500:5500]
        chain=[block(0,[(7,3500)],0,'+',len(q)),block(3500,[(7,500),(2,200),(7,2800)],2000,'+',len(q))]
        reversed_chain=[]
        for old in reversed(chain):
            body=[(op,n) for op,n in old.cigar if op not in (4,5)]
            reversed_chain.append(block(len(q)-old.qe,body,old.ref_start,'-',len(q)))
        result=detect_structures(reversed_chain,'-',{0,1},'q',reverse_complement(q),Reference(self.ref))
        parent=next(c for c in result.calls if c.type=='CPX')
        self.assertEqual((parent.info['CPX_COMPONENTS'],parent.info['NET_LEN']),('DUP,DEL',1300))

    def test_restored_branch_block_separates_two_inversions(self):
        q=self.ref[:2000]+reverse_complement(self.ref[2000:3000])+self.ref[3000:3500]+reverse_complement(self.ref[3500:4500])+self.ref[4500:6500]
        chain=[block(0,[(7,2000)],0,'+',len(q)),block(2000,[(7,1000)],2000,'-',len(q)),
            block(3000,[(7,500)],3000,'+',len(q)),block(3500,[(7,1000)],3500,'-',len(q)),block(4500,[(7,2000)],4500,'+',len(q))]
        result=detect_structures(chain,'+',{0,4},'q',q,Reference(self.ref))
        self.assertEqual([c.type for c in result.calls],['INV','INV']);self.assertEqual(result.cigar_scans,5)

    def test_short_exact_filler_separates_inv_and_del(self):
        q=self.ref[:2000]+reverse_complement(self.ref[2000:3000])+self.ref[3000:3032]+self.ref[3232:5232]
        chain=[block(0,[(7,2000)],0,'+',len(q)),block(2000,[(7,1000)],2000,'-',len(q)),block(3032,[(7,2000)],3232,'+',len(q))]
        r=detect_structures(chain,'+',{0,2},'q',q,Reference(self.ref))
        self.assertEqual(sorted(c.type for c in r.calls),['DEL','INV'])
        self.assertTrue(all(c.info.get('ROLE')!='COMPONENT' for c in r.calls))

    def test_sequence_gaps_and_N_are_not_forced_replacement(self):
        self.assertEqual(check_sequence('N'*100,self.ref[:100]).status,'ASSEMBLY_GAP')
        self.assertEqual(check_sequence(self.ref[:100],'N'*100).status,'REFERENCE_GAP')
        result=check_sequence(self.ref[1500:1600],self.ref[:3000])
        self.assertEqual(len(result.homologies),1)
        self.assertEqual([(x['query'],x['reference']) for x in result.changes],[([0,0],[0,1500]),([100,100],[1600,3000])])

    def test_delins_is_not_normalized_to_insertion(self):
        self.assertEqual(SVCandidate('DELINS','chr1',200,0,'q',{'END':400}).type,'DELINS')

    def test_junction_projection_groups_complete_pairs_and_keeps_mates(self):
        def bundle(haplotype,reverse=False):
            left=dict(chrom='chr1',position=2000,side='R')
            right=dict(chrom='chr1',position=2500,side='L')
            edge=dict(junction_id=haplotype,left=right if reverse else left,right=left if reverse else right)
            contexts=[dict(context_id=haplotype+str(i),event_id='p'+str(i),resolution='RESOLVED',junctions=[edge]) for i in range(2)]
            return dict(query='q',haplotype=haplotype,contexts=contexts,atoms=[])
        first=junction_candidates(bundle('hap1'),{})
        second=junction_candidates(bundle('hap2',True),{})
        self.assertEqual(len(first),2);self.assertEqual(first[0].info['PARENT_EVENT'],'p0,p1')
        result=cluster_and_genotype(first,second);ids={c.id for c in result}
        self.assertEqual(len(result),2)
        self.assertTrue(all(c.mateid in ids and c.genotype=='1|1' for c in result))

    def test_identical_insertions_do_not_transfer_copy_proof_between_haplotypes(self):
        first=SVCandidate('INS','chr1',2500,100,'h1',dict(SEQ='A'*100,TOPOLOGY='DUP',SOURCE_STATE='RETAINED',
            SOURCE_CHR='chr1',SOURCE_START=2400,SOURCE_END=2500,DUP_ORIENT='DIRECT',DUP_MODE='TANDEM'))
        second=SVCandidate('INS','chr1',2500,100,'h2',dict(SEQ='A'*100,SOURCE_STATE='UNKNOWN'))
        result=cluster_and_genotype([first],[second])
        self.assertEqual(len(result),1);self.assertEqual(result[0].genotype,'1|1')
        self.assertEqual(result[0].info['SOURCE_STATE'],'UNKNOWN')
        self.assertNotIn('TOPOLOGY',result[0].info)

if __name__=='__main__':unittest.main()
