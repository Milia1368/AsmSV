import random
import unittest
from asmsv_next.sequence_check import check_sequence,sparse_unique_chain,reverse_complement as rc,MAX_LINEAR_BASES


class SparseResidualTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(5677)
        self.ref=''.join(rng.choice('ACGT') for _ in range(150000))
        self.ins=''.join(rng.choice('ACGT') for _ in range(80000))

    def test_large_mixed_edits_reconstruct_complete_allele_on_both_orientations(self):
        ref=self.ref;query=ref[:1000]+self.ins+ref[1000:50000]+ref[110000:]
        for q,r in ((query,ref),(rc(query),rc(ref))):
            result=check_sequence(q,r)
            self.assertEqual(result.method,'SPARSE_EXACT_KMER');self.assertFalse(result.issues)
            self.assertEqual(len(result.changes),2)
            pieces=[];cursor=0
            for change in result.changes:
                a,b=change['reference'];s,t=change['query']
                self.assertEqual(change['type'],'GAP')
                self.assertTrue(a==b or s==t)
                pieces.extend([r[cursor:a],q[s:t]]);cursor=b
            pieces.append(r[cursor:]);self.assertEqual(''.join(pieces),q)

    def test_no_homology_repeats_and_scan_budget_abstain(self):
        self.assertTrue(check_sequence(self.ref[:60000],self.ins[:60000]).issues)
        self.assertTrue(check_sequence('ACGT'*30000,'TGCA'*30000).issues)
        self.assertEqual(sparse_unique_chain('A'*(MAX_LINEAR_BASES+1),'C'*16)[1],'SEED_BUDGET')


if __name__=='__main__':unittest.main()
