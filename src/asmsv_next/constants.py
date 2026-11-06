"""Fixed, versioned policies; these are not preset overrides.

Values have an algorithmic role. Fixing them narrows unsupported independent
tuning; it does not claim they are irrelevant on every possible input.
"""

MATCH_REWARD = 1.0
REWARD_FLOOR = 1.0
MAPQ_FULL_WEIGHT_MIN = 20
MAPQ_LOW_WEIGHT = 0.8
REFINE_MIN_ANCHOR = 16
DP_TERMINAL_RATE = 0.5
DP_EXCESS_OVERLAP_START = 10000
DP_REFERENCE_RATE = 0.1
DP_REFERENCE_CAP = 500.0
DP_WEAK_TRA_PREMIUM_RATIO = 0.2


def fixed_policy():
    return dict(match_reward=MATCH_REWARD,reward_floor=REWARD_FLOOR,
        mapq_full_weight_min=MAPQ_FULL_WEIGHT_MIN,mapq_low_weight=MAPQ_LOW_WEIGHT,
        refine_min_anchor=REFINE_MIN_ANCHOR,terminal_rate=DP_TERMINAL_RATE,
        excess_overlap_start=DP_EXCESS_OVERLAP_START,reference_rate=DP_REFERENCE_RATE,
        reference_cap=DP_REFERENCE_CAP,weak_TRA_premium_ratio=DP_WEAK_TRA_PREMIUM_RATIO,
        structural_unit='structural_gap_min',admission='known_legal_positive_reward',
        overlap='local_max',topology='additive',TRA_surcharge='information',
        aq_density='unfloored_reward_over_capped_gap_evaluation_length_v2',
        exact_linear_base_budget=2000000,sparse_seed_point_budget=20000,sparse_seed_mask=127)
