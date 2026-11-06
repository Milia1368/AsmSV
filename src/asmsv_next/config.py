class SVConfig:
    """
    Global parameter configuration for the Long-read Assembly SV detection tool.
    Centralized management of all thresholds for filtering, refinement, clustering,
    and complex SV identification.
    """
    MIN_SV_SIZE = 30
    MIN_INV_SIZE = None
    MAX_INDEL_SIZE = -1
    MAX_SV_SIZE = 1500000
    MIN_TRA_SIZE = 150000  # Legacy name: reference-gap bound for heuristic fallback, not a TRA definition.
    MIN_MAPQ = 10
    MAX_ALLELE_LEN = 10000
    ENABLE_LARGE_SEQUENCE = False
    MIN_CSV_SIZE = 400
    CLUSTER_POS_TOL = 3
    MERGE_INDEL_TOL = 3
    CLUSTER_del_LEN_TOL = 2
    CLUSTER_ins_LEN_TOL = 2
    MIN_SEQ_SIMILARITY = 1.0
    GAP_TOLERANCE = 100
    INSDEL_TOLERANCE = 10
    TANDEM_DUP_INS_THRESHOLD: int = 1000
    CALLER_PRESET = None
    DETECTOR = 'structured'
    REFERENCE_FASTA = None
    STRUCTURE_DIR = None
    SEED_MAX_BASES = 20000
    SEED_MAX_BLOCKS = 1000
    REPLACEMENT_REALIGN = 'off'
    COPY_SOURCE_SEARCH = 'adjacent'
    COPY_SOURCE_RADIUS = 100000
    VCF_INFO = 'minimal'
    GENOTYPING = 'assembly'
    SAMPLE_NAME = 'SAMPLE'
    SEQUENCE_TRANSPORT = 'auto'
    EXACT_INVERSE_MAX_BASES = 40000
    EXACT_INVERSE_MIN_COVERAGE = 0.8

    @classmethod
    def min_inv_size(cls):
        return cls.MIN_SV_SIZE if cls.MIN_INV_SIZE is None else cls.MIN_INV_SIZE
