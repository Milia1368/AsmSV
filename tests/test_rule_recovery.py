"""Counterexamples for exact structural recovery, independent of mapper splits."""
import random
import unittest

import test_structure as support
from asmsv_next.config import SVConfig
from asmsv_next.presets import resolve_preset
from asmsv_next.sequence_check import exact_inverse_edit, reverse_complement as rc
from asmsv_next.structural_detection import detect_structures


def dna(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        SVConfig.CALLER_PRESET=resolve_preset('minimap2')
        SVConfig.MIN_SV_SIZE=30;SVConfig.MIN_MAPQ=10;SVConfig.MERGE_INDEL_TOL=3
        SVConfig.COPY_SOURCE_SEARCH='adjacent';SVConfig.REPLACEMENT_REALIGN='off'

    def reverse_case(self, q, chain):
        return rc(q),[support.block(len(q)-b.qe,
            [(op,n) for op,n in b.cigar if op not in (4,5)],b.ref_start,
            '+' if b.strand=='-' else '-',len(q)) for b in reversed(chain)]

    def test_inverse_one_edit_reconstructs_both_strands(self):
        for seed in range(5):
            rng=random.Random(752+seed);ref=dna(rng,7000)
            for kind in ('DEL','INS','DELINS'):
                a,b=2500,2500+(100 if kind!='INS' else 0)
                novel=dna(rng,100) if kind!='DEL' else ''
                query=ref[:2000]+rc(ref[2000:a]+novel+ref[b:3100])+ref[3100:5100]
                chain=[support.block(0,[(7,2000)],0,'+',len(query)),
                       support.block(len(query)-2000,[(7,2000)],3100,'+',len(query))]
                for reverse in (False,True):
                    q,blocks=self.reverse_case(query,chain) if reverse else (query,chain)
                    r=detect_structures(blocks,'-' if reverse else '+',{0,1},'q',q,support.Reference(ref))
                    parent=next(c for c in r.calls if c.type=='CPX')
                    child=next(c for c in r.calls if c.type==kind)
                    self.assertEqual(parent.info['CPX_COMPONENTS'],'INV,DEL,INS' if kind=='DELINS' else 'INV,'+kind)
                    self.assertEqual(child.info['RECHECK_METHOD'],'EXACT_INVERSE')
                    s,t=parent.pos,parent.end
                    edited=ref[s:child.pos]+child.info.get('SEQ','')+ref[child.end:t]
                    self.assertEqual(ref[:s]+rc(edited)+ref[t:5100],query)
                    self.assertEqual(r.cigar_scans,2)

    def test_inverse_budget_and_negative_controls(self):
        rng=random.Random(557);ref=dna(rng,11000)
        q=rc(ref[:5000]+ref[6000:])
        self.assertIsNotNone(exact_inverse_edit(q,ref))
        self.assertIsNone(exact_inverse_edit(q,ref,max_bases=20000))
        for query,reference in [(ref,ref),(rc(ref),ref),
            (ref[:5000]+ref[6000:],ref),(q+'N',ref),
            (rc(ref[:1000]+ref[10000:]),ref),('A'*100+'C'*30+'A'*100,'T'*230)]:
            self.assertIsNone(exact_inverse_edit(query,reference))

    def test_short_separator_anchor_absorption_reconstructs_allele(self):
        for left,right in ((1,0),(0,1),(4,7),(15,15)):
            rng=random.Random(178+left+right);ref=dna(rng,7000)
            island=ref[3000:3016]
            ref=ref[:2000]+island[:left]+ref[2000+left:]
            novel=dna(rng,1000-right)+(island[-right:] if right else '')
            query=ref[:2000]+island+novel+ref[3016:5016]
            chain=[support.block(0,[(7,2000+left)],0,'+',len(query)),
                   support.block(3016-right,[(7,2000+right)],3016-right,'+',len(query))]
            for reverse in (False,True):
                q,blocks=self.reverse_case(query,chain) if reverse else (query,chain)
                r=detect_structures(blocks,'-' if reverse else '+',{0,1},'q',q,support.Reference(ref))
                self.assertEqual(sorted(c.type for c in r.calls),['DEL','INS'])
                self.assertEqual(r.cigar_scans,2)
                d=next(c for c in r.calls if c.type=='DEL');i=next(c for c in r.calls if c.type=='INS')
                reconstructed=ref[:d.pos]+ref[d.end:i.pos]+i.info['SEQ']+ref[i.pos:5016]
                self.assertEqual(reconstructed,query)

    def test_equivalent_inverted_insertion_port_preserves_original_atom(self):
        for shift in (1,2,7,15):
            rng=random.Random(982+shift);ref=dna(rng,7000);p=2000;length=500
            # Right-shift allele equivalence: inserted prefix equals target prefix.
            tail=rc(ref[p+shift:p+2*shift])
            ref=ref[:p]+tail+ref[p+shift:]
            normalized=rc(ref[p+shift:p+shift+length])
            original=ref[p:p+shift]+normalized[:-shift]
            query=ref[:p]+original+ref[p:5000]
            chain=[support.block(0,[(7,p),(1,length),(7,5000-p)],0,'+',len(query))]
            for mode in ('adjacent','local'):
              SVConfig.COPY_SOURCE_SEARCH=mode
              for reverse in (False,True):
                q,blocks=self.reverse_case(query,chain) if reverse else (query,chain)
                r=detect_structures(blocks,'-' if reverse else '+',{0},'q',q,support.Reference(ref))
                parent=next(c for c in r.calls if c.type=='DUP')
                self.assertEqual(parent.info['TARGET_POS'],p+shift)
                self.assertEqual(parent.info['DUP_NORMALIZATION'],'INS_ANCHOR_SHIFT')
                raw=next(c for c in r.calls if c.type=='INS')
                self.assertEqual((raw.pos,raw.info['SEQ']),(p,original))
                self.assertEqual(ref[:p]+raw.info['SEQ']+ref[p:5000],query)
                self.assertEqual(r.cigar_scans,1)

    def test_partial_inverse_alignment_is_completed_only_when_compatible(self):
        rng=random.Random(520);ref=dna(rng,7000);novel=dna(rng,100)
        query=ref[:2000]+rc(ref[2600:3100])+novel+rc(ref[2000:2500])+ref[3100:5100]
        for reverse in (False,True):
            chain=[support.block(0,[(7,2000)],0,'+',len(query)),
                   support.block(2600,[(7,500)],2000,'-',len(query)),
                   support.block(3100,[(7,2000)],3100,'+',len(query))]
            q,blocks=self.reverse_case(query,chain) if reverse else (query,chain)
            result=detect_structures(blocks,'-' if reverse else '+',{0,2},'q',q,support.Reference(ref))
            self.assertEqual(sorted(c.type for c in result.calls),['CPX','DELINS'])
            self.assertEqual(result.cigar_scans,3)
        chain[1]=support.block(2600,[(7,500)],2001,'-',len(query))
        result=detect_structures(chain,'+',{0,2},'q',query,support.Reference(ref))
        self.assertFalse(any(c.type=='CPX' for c in result.calls))

    def test_ambiguous_separator_does_not_force_two_independent_events(self):
        rng=random.Random(180);ref=dna(rng,7000);island=ref[3000:3016]
        ref=ref[:2000]+island[:1]+ref[2001:]
        novel=dna(rng,999)+island[-1:]
        q=ref[:2000]+island+novel+ref[3016:5016]
        chain=[support.block(0,[(7,2001)],0,'+',len(q)),
               support.block(3015,[(7,2001)],3015,'+',len(q))]
        result=detect_structures(chain,'+',{0,1},'q',q,support.Reference(ref))
        self.assertEqual([c.type for c in result.calls],['DELINS'])

    def test_source_sequence_without_retained_alignment_does_not_prove_dup(self):
        rng=random.Random(333);ref=dna(rng,15000);p=2000
        q=ref[:p]+ref[10000:10500]+ref[p:5000]
        chain=[support.block(0,[(7,p),(1,500),(7,3000)],0,'+',len(q))]
        SVConfig.COPY_SOURCE_SEARCH='local'
        result=detect_structures(chain,'+',{0},'q',q,support.Reference(ref))
        self.assertEqual([c.type for c in result.calls],['INS'])
        self.assertEqual(result.calls[0].info['SOURCE_STATE'],'UNKNOWN')

    def test_random_double_gaps_do_not_create_inverse_events(self):
        for seed in range(100):
            rng=random.Random(8010+seed)
            self.assertIsNone(exact_inverse_edit(dna(rng,1000),dna(rng,1000)))
