import unittest
from asmsv_next.alignment import AlignedSegment
from asmsv_next.config import SVConfig
from asmsv_next.presets import resolve_preset
def block(cigar,tags=None,**changes):
 data=dict(q_name='q',ref_name='chr1',ref_start=1000,cigar=cigar,mapq=60,tags=tags or {})
 data.update(changes);data['ref_end']=1000+sum(n for op,n in cigar if op in (0,2,3,7,8))
 return AlignedSegment(data)
class CigarFirstTests(unittest.TestCase):
 def setUp(self):
  SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_MAPQ=10
 def test_EQX_score_is_independent_of_aggregate_NM_nn(self):
  c=[(7,100),(1,60),(7,100)];a=block(c)
  for tags in ({'NM':63,'nn':3},{'NM':-1},{'NM':'bad'},{'NM':999999,'nn':-3}):
   with self.subTest(tags=tags):
    b=block(c,tags);self.assertTrue(b.is_valid);self.assertEqual(b.summary,a.summary)
    self.assertEqual(dict(b.score_root.raw.tags),tags)
 def test_complete_EQX_does_not_parse_redundant_difference_tags(self):
  a=block([(7,100),(8,10),(7,100)],{'MD':'malformed','cs':':210','NM':0})
  self.assertTrue(a.is_valid);self.assertEqual((a.summary.matches,a.summary.mismatches),(200,10))
  self.assertEqual(a.as_score,10)
 def test_inserted_unknowns_cannot_cancel_gap_cost(self):
  c=[(7,100),(1,5000),(7,100)]
  a=block(c,{'NM':5000,'nn':5000})
  self.assertEqual(a.summary.matches,200)
  self.assertAlmostEqual(a.as_score,200-SVConfig.CALLER_PRESET.aq.gap_cost(5000))
 def test_M_requires_difference_evidence(self):
  a=block([(0,200)]);self.assertFalse(a.is_valid);self.assertIsNone(a.summary.matches)
 def test_M_MD_is_positional_and_does_not_need_NM(self):
  a=block([(0,100),(1,60),(0,100)],{'MD':'100A99'})
  self.assertTrue(a.is_valid);self.assertTrue(a.summary.positional)
  self.assertEqual((a.summary.matches,a.summary.mismatches),(199,1))
 def test_M_MD_has_priority_over_inconsistent_NM(self):
  for nm in (-1,63,1000,'bad'):
   a=block([(0,100),(1,60),(0,100)],{'MD':'100A99','NM':nm})
   self.assertTrue(a.is_valid);self.assertEqual((a.summary.matches,a.summary.mismatches),(199,1))
 def test_M_cs_has_priority_over_inconsistent_NM(self):
  a=block([(0,100),(1,60),(0,100)],{'cs':':100+'+'a'*60+'*ac:99','NM':1000})
  self.assertTrue(a.is_valid);self.assertEqual((a.summary.matches,a.summary.mismatches),(199,1))
 def test_conflicting_M_positional_tags_still_reject(self):
  a=block([(0,200)],{'MD':'100A99','cs':':200','NM':1})
  self.assertFalse(a.is_valid);self.assertIn('MD_CS_CONFLICT',a.rejection_reasons)
 def test_malformed_M_positional_tag_still_rejects(self):
  a=block([(0,200)],{'MD':'bad','NM':1})
  self.assertFalse(a.is_valid);self.assertIn('INVALID_MD',a.rejection_reasons)
 def test_M_NM_fallback_subtracts_explicit_X_and_I_D_once(self):
  a=block([(0,100),(8,2),(1,50),(2,20),(0,100)],{'NM':75})
  self.assertTrue(a.is_valid);self.assertFalse(a.summary.positional)
  self.assertEqual((a.summary.matches,a.summary.mismatches),(197,5))
 def test_M_NM_only_does_not_invent_exact_slice_positions(self):
  a=block([(0,200)],{'NM':2});b=a.split_at_query_coords(30,170)
  self.assertTrue(a.is_valid);self.assertFalse(b.is_valid)
  self.assertIn('UNLOCATED_SLICE_DIFFERENCES',b.rejection_reasons)
 def test_invalid_M_NM_is_rejected_when_it_is_required(self):
  for nm in (-1,201,'3'):
   self.assertFalse(block([(0,200)],{'NM':nm}).is_valid)
 def test_refskip_and_MAPQ_gates_remain(self):
  self.assertFalse(block([(7,100),(3,50),(7,100)]).is_valid)
  self.assertFalse(block([(7,200)],mapq=255).is_valid)
 def test_EQX_and_lossless_M_MD_have_identical_scores_and_slices(self):
  for reverse in (False,True):
   for insertion in (1,49,50,51,500):
    a=block([(7,40),(8,2),(7,58),(1,insertion),(7,100)],is_reverse=reverse)
    b=block([(0,100),(1,insertion),(0,100)],{'MD':'40A0C158','NM':insertion+9,'nn':7},is_reverse=reverse)
    self.assertTrue(a.is_valid and b.is_valid);self.assertEqual(a.summary,b.summary)
    self.assertEqual(a.as_score,b.as_score)
    for left,right in ((10,80),(20,180+insertion)):
     x=a.split_at_query_coords(left,right);y=b.split_at_query_coords(left,right)
     self.assertEqual(x.summary,y.summary);self.assertEqual((x.qs,x.qe,x.ref_start,x.ref_end),(y.qs,y.qe,y.ref_start,y.ref_end))
if __name__=='__main__':unittest.main()
