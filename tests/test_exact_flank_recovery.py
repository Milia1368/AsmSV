import random
import unittest
from asmsv_next.config import SVConfig
from asmsv_next.presets import resolve_preset
from asmsv_next.sequence_check import exact_flank_model,check_sequence,reverse_complement as rc
from asmsv_next.structural_detection import detect_structures
from test_structure import Reference,block


class ExactFlankRecoveryTests(unittest.TestCase):
    def setUp(self):
        SVConfig.CALLER_PRESET=resolve_preset('minimap2')
        SVConfig.MIN_SV_SIZE=50;SVConfig.MIN_MAPQ=10
        SVConfig.COPY_SOURCE_SEARCH='off';SVConfig.REPLACEMENT_REALIGN='off'
        rng=random.Random(270611)
        self.ref=''.join(rng.choice('ACGT') for _ in range(165000))

    def test_large_flanked_inverse_bypasses_seed_budget_with_exact_reconstruction(self):
        r=self.ref[:100]+self.ref[100:100100]+self.ref[100100:100200]
        q=r[:100]+rc(r[100:100100])+r[100100:]
        result=check_sequence(q,r,max_bases=20000)
        self.assertFalse(result.issues);self.assertEqual(len(result.changes),1)
        c=result.changes[0];s,e=c['reference']
        self.assertEqual(c['type'],'INV')
        self.assertEqual(r[:s]+rc(r[s:e])+r[e:],q)

    def test_short_and_long_backbone_inversions_and_global_reverse_complement(self):
        for length in (50,54,59,99,500,33891,100000):
            query=self.ref[:1000]+rc(self.ref[1000:1000+length])+self.ref[1000+length:2000+length]
            chain=[block(0,[(7,1000)],0,'+',len(query)),
                   block(1000+length,[(7,1000)],1000+length,'+',len(query))]
            for reverse in (False,True):
                if reverse:
                    transformed=[block(len(query)-b.qe,[(7,b.qe-b.qs)],b.ref_start,'-',len(query)) for b in reversed(chain)]
                    result=detect_structures(transformed,'-',{0,1},'q',rc(query),Reference(self.ref))
                else:result=detect_structures(chain,'+',{0,1},'q',query,Reference(self.ref))
                calls=[c for c in result.calls if c.type=='INV']
                self.assertEqual(len(calls),1,(length,reverse))
                c=calls[0]
                self.assertEqual(self.ref[:c.pos]+rc(self.ref[c.pos:c.end])+self.ref[c.end:2000+length],query)
                self.assertGreaterEqual(c.sv_len,50)

    def test_boundary_absorbed_54bp_inverse_recovers_from_a_48bp_gap(self):
        # Three bases at both true inversion boundaries equal their counterparts.
        middle=self.ref[3000:3048];core='ACG'+middle+'CGT'
        ref=self.ref[:1000]+core+self.ref[1054:2054]
        query=ref[:1000]+rc(core)+ref[1054:]
        chain=[block(0,[(7,1003)],0,'+',len(query)),
               block(1051,[(7,1003)],1051,'+',len(query))]
        result=detect_structures(chain,'+',{0,1},'q',query,Reference(ref))
        inv=next(c for c in result.calls if c.type=='INV')
        self.assertGreaterEqual(inv.sv_len,50)
        self.assertEqual(ref[:inv.pos]+rc(ref[inv.pos:inv.end])+ref[inv.end:],query)

    def test_below_minimum_and_palindromic_no_change_are_not_inversions(self):
        for length in (16,25):
            query=self.ref[:1000]+rc(self.ref[1000:1000+length])+self.ref[1000+length:2000+length]
            chain=[block(0,[(7,1000)],0,'+',len(query)),
                   block(1000+length,[(7,1000)],1000+length,'+',len(query))]
            result=detect_structures(chain,'+',{0,1},'q',query,Reference(self.ref))
            self.assertFalse(any(c.type=='INV' for c in result.calls))
        self.assertEqual(exact_flank_model('ACGT','ACGT').changes,[])
        self.assertIsNone(exact_flank_model('N'*100,'N'*100))

    def test_large_exact_indel_flanks_restore_only_the_uncovered_edit(self):
        ref=self.ref[:50000]
        for query,expected in ((ref[:100]+'T'*40000+ref[100:],'GAP'),
                               (ref[:100]+ref[40100:],'GAP')):
            result=check_sequence(query,ref,max_bases=20000)
            self.assertFalse(result.issues)
            self.assertEqual([c['type'] for c in result.changes],[expected])
            c=result.changes[0];qs,qe=c['query'];rs,re=c['reference']
            self.assertEqual(ref[:rs]+query[qs:qe]+ref[re:],query)

    def test_unexplained_or_N_containing_large_gap_is_not_forced_to_a_model(self):
        self.assertIsNone(exact_flank_model(self.ref[:1000],self.ref[3000:4000]))
        self.assertEqual(check_sequence('N'+self.ref[:30000],self.ref[:30001]).status,'ASSEMBLY_GAP')


if __name__=='__main__':unittest.main()
