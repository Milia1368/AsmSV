import contextlib
import io
import random
import unittest
from asmsv_next.cli import arguments
from asmsv_next.config import SVConfig
from asmsv_next.structures import SVCandidate,apply_indel_size_policy
from asmsv_next.structural_detection import detect_structures
from asmsv_next.sequence_check import reverse_complement as rc
from asmsv_next.presets import resolve_preset
from test_structure import Reference,block
class SizePolicyTests(unittest.TestCase):
 def setUp(self):
  keys=('MIN_SV_SIZE','MIN_INV_SIZE','MAX_INDEL_SIZE','MAX_ALLELE_LEN','COPY_SOURCE_SEARCH','CALLER_PRESET','MIN_MAPQ')
  original={k:getattr(SVConfig,k) for k in keys}
  self.addCleanup(lambda:[setattr(SVConfig,k,v) for k,v in original.items()])
  SVConfig.MIN_SV_SIZE=50;SVConfig.MIN_INV_SIZE=None;SVConfig.MAX_INDEL_SIZE=-1
  SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.COPY_SOURCE_SEARCH='off';SVConfig.MIN_MAPQ=10
 def inversion(self,length,split=False,reverse=False):
  rng=random.Random(834+length);core='A'+''.join(rng.choices('ACGT',k=length-2))+'A'
  ref='C'*100+core+'C'*100;query=ref[:100]+rc(core)+ref[100+length:]
  chain=[block(0,[(7,100)],0,'+',len(query)),block(100+length,[(7,100)],100+length,'+',len(query))]
  if split:chain.insert(1,block(100,[(7,length)],100,'-',len(query)))
  if reverse:
   chain=[block(len(query)-b.qe,[(7,b.qe-b.qs)],b.ref_start,'-' if b.strand=='+' else '+',len(query)) for b in reversed(chain)]
   query=rc(query)
  return detect_structures(chain,'-' if reverse else '+',{0,len(chain)-1},'q',query,Reference(ref))
 def test_default_min_inv_follows_min_sv_at_inclusive_boundary(self):
  for n in (49,50,51):
   result=self.inversion(n)
   self.assertEqual(any(c.type=='INV' for c in result.calls),n>=50)
 def test_min_inv_is_independent_for_gap_and_split_evidence_on_both_strands(self):
  for split in (False,True):
   for reverse in (False,True):
    SVConfig.MIN_INV_SIZE=30
    result=self.inversion(40,split,reverse)
    inv=[c for c in result.calls if c.type=='INV']
    self.assertEqual(len(inv),1);self.assertEqual(inv[0].sv_len,40)
    SVConfig.MIN_INV_SIZE=41
    result=self.inversion(40,split,reverse)
    self.assertFalse(any(c.type=='INV' for c in result.calls))
    self.assertFalse(any(c.get('junctions') for c in result.contexts if c['frame']=='BELOW_MIN_INV'))
 def test_optional_max_keeps_types_and_is_inclusive(self):
  SVConfig.MAX_INDEL_SIZE=10000
  calls=[SVCandidate(k,'chr1',10,n,'q',{'SEQ':'A'*n} if k=='INS' else {}) for k in ('INS','DEL') for n in (9999,10000,10001)]
  apply_indel_size_policy(calls)
  self.assertEqual([c.filter for c in calls],['PASS','PASS','OVERSIZE']*2)
  self.assertEqual([c.type for c in calls],['INS']*3+['DEL']*3)
 def test_standalone_INV_policy_does_not_disable_CPX_inversion_frame(self):
  rng=random.Random(998);ref=''.join(rng.choices('ACGT',k=12000))
  query=ref[:2000]+rc(ref[2700:3500])+rc(ref[2000:2500])+ref[3500:5500]
  chain=[block(0,[(7,2000)],0,'+',len(query)),
   block(2000,[(7,500),(2,200),(7,800)],2000,'-',len(query)),
   block(3300,[(7,2000)],3500,'+',len(query))]
  SVConfig.MIN_INV_SIZE=2000
  result=detect_structures(chain,'+',{0,2},'q',query,Reference(ref))
  parents=[c for c in result.calls if c.type=='CPX']
  self.assertEqual(len(parents),1);self.assertEqual(parents[0].info['FRAME'],'INV')
 def test_unlimited_policy_and_unrelated_types(self):
  calls=[SVCandidate(k,'chr1',10,200000,'q') for k in ('INS','DEL','INV','DUP','BND')]
  apply_indel_size_policy(calls);self.assertTrue(all(c.filter=='PASS' for c in calls))
  SVConfig.MAX_INDEL_SIZE=10000;apply_indel_size_policy(calls)
  self.assertEqual([c.filter for c in calls],['OVERSIZE','OVERSIZE','PASS','PASS','PASS'])
 def test_optional_max_preserves_existing_filter(self):
  SVConfig.MAX_INDEL_SIZE=10000
  c=SVCandidate('INS','chr1',10,10001,'q',filter_str='PARTIAL')
  apply_indel_size_policy([c]);apply_indel_size_policy([c])
  self.assertEqual(c.filter,'OVERSIZE;PARTIAL')
 def test_INS_ALT_is_not_truncated_at_any_display_or_reporting_limit(self):
  SVConfig.MAX_ALLELE_LEN=1;SVConfig.MAX_INDEL_SIZE=10000
  for n in (9999,10000,10001,250000):
   sequence=('ACGT'*((n+3)//4))[:n]
   c=SVCandidate('INS','chr1',10,n,'q',{'SEQ':sequence})
   apply_indel_size_policy([c]);record=c.to_vcf_record(1,reference=Reference('C'*100))
   self.assertEqual(record[4],'C'+sequence);self.assertEqual(c.sv_len,n)
 def test_cli_defaults_aliases_and_invalid_bounds(self):
  base=['--aligner','minimap2','--bam1','h1.bam','--bam2','h2.bam','--ref','r.fa']
  a=arguments(base);self.assertEqual((a.min_inv,a.max_indel_size,a.max_ref_gap),(None,-1,150000))
  a=arguments(base+['--min-inv-size','40','--max-indel-size','10000','--max-ref-gap','500000'])
  self.assertEqual((a.min_inv,a.max_indel_size,a.max_ref_gap),(40,10000,500000))
  for extra in (['--min-inv','0'],['--max-indel-size','0'],['--max-ref-gap','-1']):
   with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):arguments(base+extra)
if __name__=='__main__':unittest.main()
