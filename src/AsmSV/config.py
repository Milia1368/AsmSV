# ==============================================================================
# Global Configuration and Default Parameters
# ==============================================================================

class SVConfig:
    """
    Global parameter configuration for the Long-read Assembly SV detection tool.
    Centralized management of all thresholds for filtering, refinement, clustering, 
    and complex SV identification.
    """
    
    # ================= Base Constants =================
    # These are "base" values. It is not recommended to modify them directly via 
    # the command line, but rather shift them overall via CLUSTER_LEN_TOL.
    _BASE_INS_LEN_TOL = 0
    _BASE_DEL_LEN_TOL = 2
    # ================= Basic Filtering Parameters =================
    MIN_SV_SIZE = 30                # Minimum SV size to report (bp)
    MAX_SV_SIZE = 1500000           # Maximum SV size (larger events may be assembly errors)
    MIN_TRA_SIZE = 150000           # Minimum gap size to classify as a Translocation (TRA/BND)
    MIN_MAPQ = 10                   # Minimum mapping quality score
    MAX_ALLELE_LEN = 5000           # Maximum sequence length to write in the ALT column
    ENABLE_ALLELE = True            # Enable writing explicit allele sequences (instead of symbolic)
    MIN_CSV_SIZE=400                # Minimum SV size to report ins-del rpl event (bp)
    MAX_CSV_SIZE=40000              # Maximum SV size to report ins-del rpl event (bp)
   
    # ================= Clustering & Merging Parameters =================
    CLUSTER_POS_TOL = 3             # Breakpoint position distance tolerance for clustering
    MERGE_INDEL_TOL = 3             # Max distance tolerance to merge adjacent identical INDELs
    # Core variable: General length tolerance offset
    CLUSTER_LEN_TOL = 0 
    CLUSTER_del_LEN_TOL = _BASE_DEL_LEN_TOL + CLUSTER_LEN_TOL
    CLUSTER_ins_LEN_TOL = _BASE_INS_LEN_TOL + CLUSTER_LEN_TOL
    MIN_SEQ_SIMILARITY = 1.0        # Minimum sequence similarity required for merging INS

    # ================= Complex SV & Gap Processing =================
    GAP_TOLERANCE = 100             # Coordinate tolerance for identifying complex structures / Tandem Duplications
    INSDEL_TOLERANCE = 10           # Length difference tolerance for calling split-read INS/DELs
    
    # ================= Density Filtering & Miscellaneous =================
    DENSITY_WINDOW = 1000           # Window size for high-density noise filtering (bp)
    MAX_SV_IN_WINDOW = 3            # Maximum number of SVs allowed within a density window
    TRA_RATIO_THRESHOLD = 5.0       # Distance ratio threshold for evaluating translocations/inversions
    ENABLE_DENSITY_FILTER = False   # Enable high-density noise filtering for false positives
    ENABLE_SPLIT_SEGMENTS = False   # Enable segment splitting/refinement before DP chaining

    # ================= Tandem / Dispersed Duplication Classification =================
    # 真正的 tandem dup 拷贝紧贴或重叠源位点（通常 ≤ 500 bp）
    # 原 condense_orphan_blocks 用 max_ref_dist=100000，分类过于宽泛
    TANDEM_MAX_REF_DIST: int = 500

    # 近端散在重复（PROXIMAL）：距 gap 500 bp~10 kb，机制不确定，单独标注
    PROXIMAL_MAX_REF_DIST: int = 10_000

    # DUP:TANDEM 的 INS/DUP SVTYPE 长度分界（参考 Manta 惯例）
    # 短拷贝（< 1 kb）→ INS；长拷贝（≥ 1 kb）→ DUP
    TANDEM_DUP_INS_THRESHOLD: int = 1_000

    # pair_inv_signals 中 InvDUP:TANDEM 判定所需的最小 backbone 倒退量
    # 原用 MIN_SV_SIZE（30-50 bp），噪声引起的误判；改为 200 bp
    MIN_INVDUP_TANDEM_BACK: int = 200
    #----backbone---
    BB_OVERLAP_MIN_BASE=200
    BB_OVERLAP_MAX_CEILING=5000
    BB_OVERLAP_RATIO_THRES=0.5
