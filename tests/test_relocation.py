import random
import unittest
from asmsv_next.relocation import link_relocations,rotation_key
from asmsv_next.structures import SVCandidate
from asmsv_next.sequence_check import reverse_complement as rc
from asmsv_next.output_views import junction_candidates
from test_structure import Reference


class RelocationTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(26788);self.ref=''.join(rng.choice('ACGT') for _ in range(10000))

    def calls(self,inverse=False):
        sequence=self.ref[1000:1800]
        deletion=SVCandidate('DEL','chr1',1000,800,'donor',{'ATOM_ID':'del'})
        insertion=SVCandidate('INS','chr2',3000,800,'target',{'SEQ':rc(sequence) if inverse else sequence,'ATOM_ID':'ins'})
        return deletion,insertion

    def test_direct_and_inverse_cross_contig_links_preserve_edits(self):
        for inverse in (False,True):
            d,i=self.calls(inverse);original=(i.pos,i.info['SEQ'],d.pos,d.end)
            relations=link_relocations([d,i],'hap1',Reference(self.ref))
            self.assertEqual(len(relations),1)
            self.assertEqual(relations[0]['source'][-1],'INVERTED' if inverse else 'DIRECT')
            self.assertEqual((i.pos,i.info['SEQ'],d.pos,d.end),original)
            self.assertEqual(i.info['SOURCE_STATE'],'REMOVED')
            bundle=dict(query='target',haplotype='hap1',contexts=[],atoms=[],relocations=relations)
            joins=junction_candidates(bundle,{})
            self.assertEqual(len(joins),4)
            self.assertEqual({c.info['RELOCATION_ID'] for c in joins},{relations[0]['relation_id']})

    def test_no_deleted_source_retained_source_or_ambiguous_sources_abstain(self):
        d,i=self.calls()
        self.assertFalse(link_relocations([i],'hap1',Reference(self.ref)))
        i.info['SOURCE_STATE']='RETAINED'
        self.assertFalse(link_relocations([d,i],'hap1',Reference(self.ref)))
        d,i=self.calls();d2=SVCandidate('DEL','chr3',1000,800,'another')
        self.assertFalse(link_relocations([d,i,d2],'hap1',Reference(self.ref)))
        self.assertEqual(i.info['SOURCE_ISSUE'],'OBSERVED_DELETION_SOURCE_AMBIGUOUS')

    def test_rotated_index_requires_actual_target_allele_equivalence(self):
        d,i=self.calls();source=i.info['SEQ'];i.info['SEQ']=source[-2:]+source[:-2]
        different=lambda base:next(c for c in 'ACGT' if c!=base)
        reference=self.ref[:999]+different(source[-1])+self.ref[1000:1800]+different(source[0])+self.ref[1801:]
        self.assertEqual(rotation_key(source),rotation_key(i.info['SEQ']))
        self.assertFalse(link_relocations([d,i],'hap1',Reference(reference)))
        # Positive shift two restores source DNA only if those two reference
        # target bases are the same prefix of the actual inserted allele.
        reference=reference[:3000]+source[-2:]+reference[3002:]
        links=link_relocations([d,i],'hap1',Reference(reference))
        self.assertEqual(len(links),1);self.assertEqual(links[0]['target_boundary_shift'],2)

    def test_donor_boundary_equivalence_is_independently_verified(self):
        d,i=self.calls();source=i.info['SEQ']
        reference=self.ref[:998]+source[-2:]+self.ref[1000:]
        i.info['SEQ']=source[-2:]+source[:-2]
        links=link_relocations([d,i],'hap1',Reference(reference))
        self.assertEqual(len(links),1)
        relation=links[0];s,e=relation['source'][1:3]
        self.assertEqual(reference[:d.pos]+reference[d.end:],reference[:s]+reference[e:])
        self.assertEqual(reference[s:e],i.info['SEQ'])


if __name__=='__main__':unittest.main()
