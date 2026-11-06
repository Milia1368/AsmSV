import unittest
from asmsv_next.presets import resolve_preset
from asmsv_next.config import SVConfig
from test_core import block


class AQDensityTests(unittest.TestCase):
    def setUp(self):SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_MAPQ=10

    def test_density_has_no_gap_threshold_jump_and_reward_is_unchanged(self):
        p=SVConfig.CALLER_PRESET.aq
        for op in (1,2):
            values=[]
            for n in (48,49,50,51,1000):
                a=block([(7,100),(op,n),(7,100)])
                reward=200-p.gap_cost(n)
                self.assertAlmostEqual(a.summary.node_reward,reward)
                self.assertAlmostEqual(a.summary.evaluation_length,200+min(n,p.structural_gap_min))
                self.assertAlmostEqual(a.aq,reward/(200+min(n,p.structural_gap_min)))
                self.assertAlmostEqual(a.as_score,max(1,reward))
                values.append(a.aq)
            self.assertTrue(all(a>b for a,b in zip(values,values[1:])))
            self.assertLess(abs(values[2]-values[1]),.006)

    def test_fractional_insertion_slice_length_and_mapq_once(self):
        a=block([(7,100),(1,1000),(7,100)],mapq=10)
        # Direct root scores retain exact proportional structural costs and
        # capped evaluation length even when a view ends within an insertion.
        summary=a.score_root.index.score(0,433)
        expected_reward=100-.333*SVConfig.CALLER_PRESET.aq.gap_cost(1000)
        self.assertAlmostEqual(summary.node_reward,expected_reward)
        self.assertAlmostEqual(summary.evaluation_length,116.65)
        self.assertAlmostEqual(summary.aq,expected_reward/116.65)
        self.assertAlmostEqual(a.aq,.8*a.aq_sequence)


if __name__=='__main__':unittest.main()
