import hashlib
from pathlib import Path
import random
import tempfile
import unittest
import pysam
from asmsv_next.allele_normalization import canonical_tandem_indels,shifted_insertion,canonical_inversion
from asmsv_next.sequence_output import SequenceStore
from asmsv_next.sequence_check import reverse_complement as rc
from asmsv_next.structures import SVCandidate
from asmsv_next.structural_detection import detect_structures
from asmsv_next.config import SVConfig
from asmsv_next.presets import resolve_preset
from test_structure import Reference,block
from asmsv_next.output_views import junction_candidates


class AllelePayloadTests(unittest.TestCase):
    def setUp(self):
        SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_SV_SIZE=50
        SVConfig.MIN_MAPQ=10;SVConfig.COPY_SOURCE_SEARCH='off';SVConfig.REPLACEMENT_REALIGN='off'
        rng=random.Random(110628)
        self.ref=''.join(rng.choice('ACGT') for _ in range(40000))

    def test_long_literal_ALT_and_deduplicated_FASTA_roundtrip(self):
        sequence=self.ref[:25000]
        a=SVCandidate('INS','chr1',100,25000,'q',{'SEQ':sequence})
        b=SVCandidate('INS','chr1',200,25000,'q',{'SEQ':sequence})
        with tempfile.TemporaryDirectory() as folder:
            store=SequenceStore(Path(folder)/'insertions.fa')
            store.annotate(a);store.annotate(b);store.finish();store.abort()
            self.assertEqual(len(store.seen),1)
            with pysam.FastaFile(str(store.path)) as fasta:
                recovered=fasta.fetch(a.info['SEQID'])
            self.assertEqual(recovered,sequence)
            self.assertEqual(hashlib.sha256(recovered.encode()).hexdigest(),a.info['SEQSHA256'])
            record=a.to_vcf_record(1,reference=Reference(self.ref))
            self.assertEqual(record[4][1:],sequence)

    def test_tandem_projection_is_exact_and_keeps_original_atom(self):
        s,t=1000,25000;payload=self.ref[s:t]
        atom=SVCandidate('INS','chr1',s,t-s,'q',{'SEQ':payload,'PARENT_EVENT':'parent'})
        parent=SVCandidate('DUP','chr1',s,t-s,'q',{'DUP_MODE':'TANDEM','TARGET_CHR':'chr1','TARGET_POS':t})
        result=canonical_tandem_indels([atom],{'parent':parent},Reference(self.ref))[0]
        self.assertEqual(atom.pos,s);self.assertIsNot(result,atom);self.assertEqual(result.pos,t)
        self.assertEqual(self.ref[:s]+payload+self.ref[s:],self.ref[:t]+result.info['SEQ']+self.ref[t:])
        self.assertIsNone(shifted_insertion(Reference(self.ref),'chr1',s,'A'*100,t))

    def test_padded_microoverlap_returns_actual_full_allele_on_both_strands(self):
        payload=self.ref[10000:35000];ref=self.ref[:3000];query=ref[:1000]+payload+ref[1000:]
        # Reference boundaries overlap by seven bases: padded comparison must
        # restore the duplicated target bases before extracting the payload.
        chain=[block(0,[(7,1000)],0,'+',len(query)),block(26000,[(7,2000)],993,'+',len(query))]
        # Make the second anchor truthful and include the overlapped seven bases
        # in the inserted allele rather than inventing an unverified alignment.
        query=ref[:1000]+payload+ref[993:]
        for reverse in (False,True):
            if reverse:
                nodes=[block(len(query)-b.qe,[(7,b.qe-b.qs)],b.ref_start,'-',len(query)) for b in reversed(chain)]
                result=detect_structures(nodes,'-',{0,1},'q',rc(query),Reference(ref))
            else:result=detect_structures(chain,'+',{0,1},'q',query,Reference(ref))
            ins=next(c for c in result.calls if c.type=='INS')
            self.assertEqual(ref[:ins.pos]+ins.info['SEQ']+ref[ins.pos:],query)

    def test_inversion_complementary_flanks_are_equivalent_and_query_verified(self):
        middle=self.ref[4000:5000];ref=self.ref[:1000]+'ACG'+middle+'CGT'+self.ref[2006:3006]
        start,end=1003,2003;query=ref[:start]+rc(ref[start:end])+ref[end:]
        normalized=canonical_inversion(Reference(ref),'chr1',start,end,query,(start,end),'+')
        self.assertIsNotNone(normalized)
        s,e,_=normalized
        self.assertEqual(ref[:s]+rc(ref[s:e])+ref[e:],query)
        self.assertIsNone(canonical_inversion(Reference(ref),'chr1',start,end,'N'*len(query),(start,end),'+'))

    def test_known_source_insertion_keeps_two_source_joins_without_copy_claim(self):
        ref=self.ref;query=ref[:1000]+ref[20000:20500]+ref[1000:2000]
        chain=[block(0,[(7,1000)],0,'+',len(query)),
               block(1000,[(7,500)],20000,'+',len(query)),
               block(1500,[(7,1000)],1000,'+',len(query))]
        result=detect_structures(chain,'+',{0,2},'q',query,Reference(ref),'hap1')
        ins=next(c for c in result.calls if c.type=='INS')
        self.assertEqual(ins.info['SOURCE_STATE'],'UNKNOWN')
        bundle=dict(query='q',haplotype='hap1',contexts=result.contexts,atoms=result.atoms)
        calls=junction_candidates(bundle,{})
        self.assertEqual(len(calls),4)
        self.assertEqual(len({c.info['JUNCTION_ID'] for c in calls}),2)
        self.assertTrue(all(c.filter=='PASS' for c in calls))


if __name__=='__main__':unittest.main()
