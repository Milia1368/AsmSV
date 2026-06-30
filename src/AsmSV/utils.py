import os
import time
import logging
import argparse
import pysam
import edlib
from typing import List, Dict, Tuple, Any
from collections import defaultdict, deque

# Import core structures
from AsmSV.structures import SVCandidate 
from AsmSV.config import SVConfig

logger = logging.getLogger(__name__)

# ==============================================================================
# Global Constants
# ==============================================================================

# Pre-compiled mapping table for reverse complement operations
COMPLEMENT_MAP = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')

# ==============================================================================
# Argument Parsing and Configuration
# ==============================================================================

def parseArgs(argv: List[str]) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Advanced Long-read Assembly SV Detector (v1.2.1)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # === Group 1: Input/Output (Required) ===
    io_group = parser.add_argument_group("Input/Output Arguments")
    io_group.add_argument("--bam1", required=True, metavar="FILE", help="Input BAM for Hap1 (sorted & indexed)")
    io_group.add_argument("--bam2", required=True, metavar="FILE", help="Input BAM for Hap2 (sorted & indexed)")
    io_group.add_argument("--ref", required=True, metavar="FILE", help="Reference Genome FASTA (indexed)")
    io_group.add_argument("--output", "-o", metavar="PATH", default="./results.vcf",
                          help="Output VCF path. If a directory is given, defaults to results.vcf inside it. [Default: ./results.vcf]")
    io_group.add_argument("--mask", metavar="FILE", default=None,
                          help="BED file of regions to exclude from output. [Default: None]")
    io_group.add_argument("--sorted-bam1", metavar="FILE", default=None,
                          help="Pre-sorted (queryname order) BAM for Hap1. Skips the sort step entirely.\n"
                               "When provided the file is never deleted regardless of --no-retain-sorted.")
    io_group.add_argument("--sorted-bam2", metavar="FILE", default=None,
                          help="Pre-sorted (queryname order) BAM for Hap2. Skips the sort step entirely.\n"
                               "When provided the file is never deleted regardless of --no-retain-sorted.")
    
    # === Group 2: VCF Format Options ===
    fmt_group = parser.add_argument_group("VCF Format Options")
    fmt_group.add_argument(
                        "--enable-large-sequence",
                        dest="enable_large_sequence",
                        action="store_true",
                        help=f"""Always write full inserted/deleted sequences in the VCF ALT/REF column,
regardless of length. By default (off), sequences longer than the
internal threshold ({SVConfig.MAX_ALLELE_LEN} bp) are represented as symbolic alleles
(e.g. <INS>, <DEL>) to keep the VCF compact.
Enable this flag when downstream tools require explicit sequences.
WARNING: may produce very large VCF files for large SVs.
NOTE: DUP/INV/CSV always use symbolic alleles regardless of this flag."""
                        )
    parser.set_defaults(enable_large_sequence=False)

    # === Group 3: Variant Filtering Parameters ===
    filter_group = parser.add_argument_group("Filter Parameters")
    filter_group.add_argument("--min-sv", type=int, default=SVConfig.MIN_SV_SIZE, metavar="INT",
                              help=f"Minimum SV span (bp) to report. Variants smaller than this are\n"
                                   f"skipped during candidate extraction (CIGAR and split-read). [Default: {SVConfig.MIN_SV_SIZE}]")
    filter_group.add_argument("--max-sv", type=int, default=SVConfig.MAX_SV_SIZE, metavar="INT",
                              help=f"Maximum SV span (bp) to report. Span SVs (DEL/INS/INV/DUP/CSV)\n"
                                   f"exceeding this size are tagged FILTER=OVERSIZE in the VCF instead\n"
                                   f"of being silently dropped, so they can be filtered post-hoc.\n"
                                   f"Set to -1 to disable the upper bound entirely. [Default: {SVConfig.MAX_SV_SIZE}]\n"
                                   f"NOTE: does not affect VCF sequence output (see --enable-large-sequence).")
    filter_group.add_argument("--min-mapq", type=int, default=SVConfig.MIN_MAPQ, metavar="INT",
                              help=f"Minimum mapping quality for an alignment to be used in SV calling. [Default: {SVConfig.MIN_MAPQ}]")
    filter_group.add_argument("--insdel-tolerance", type=int, default=SVConfig.INSDEL_TOLERANCE, metavar="INT",
                              help=f"Length difference tolerance (bp) between query-gap and ref-gap when\n"
                                   f"classifying a split-read signal as INS vs DEL. [Default: {SVConfig.INSDEL_TOLERANCE}]")
    filter_group.add_argument("--min-tra-size", type=int, default=SVConfig.MIN_TRA_SIZE, metavar="INT",
                              help=f"Minimum inter-alignment gap (bp) on the reference to classify an event\n"
                                   f"as a translocation (TRA/BND) instead of a local SV. [Default: {SVConfig.MIN_TRA_SIZE}]")
    filter_group.add_argument("--min-csv", type=int, default=SVConfig.MIN_CSV_SIZE, metavar="INT",
                              help=f"Minimum size (bp) for a combined INS+DEL replacement event (CSV/INS_DEL)\n"
                                   f"to be reported. [Default: {SVConfig.MIN_CSV_SIZE}]")
    filter_group.add_argument("--max-csv", type=int, default=SVConfig.MAX_CSV_SIZE, metavar="INT",
                              help=f"Maximum size (bp) for a CSV/INS_DEL replacement event. Events exceeding\n"
                                   f"this limit are skipped during extraction. [Default: {SVConfig.MAX_CSV_SIZE}]")
    filter_group.add_argument("--split-segments", action="store_true",
                              help="Enable query-interval overlap resolution before DP chaining.\n"
                                   "Splits overlapping alignment segments to improve SV boundary accuracy.\n"
                                   "[Default: off]")
    
    # === Group 4: Clustering Parameters ===
    cluster_group = parser.add_argument_group("Clustering Parameters")
    cluster_group.add_argument("--merge-indels-tol", type=int, default=SVConfig.MERGE_INDEL_TOL, metavar="INT",
                               help=f"Distance tolerance to merge adjacent identical variants. [Default: {SVConfig.MERGE_INDEL_TOL}]")
    cluster_group.add_argument("--cluster-pos-tol", type=int, default=SVConfig.CLUSTER_POS_TOL, metavar="BP",
                               help=f"Max breakpoint distance tolerance for merging. [Default: {SVConfig.CLUSTER_POS_TOL}]")
    cluster_group.add_argument("--cluster-len-tol", type=int, default=SVConfig.CLUSTER_LEN_TOL, metavar="BP",
                               help=f"Length tolerance adjustment. [Default: {SVConfig.CLUSTER_LEN_TOL}]")
    cluster_group.add_argument("--min-seq-sim", type=float, default=SVConfig.MIN_SEQ_SIMILARITY, metavar="FLOAT",
                               help=f"Min sequence similarity for INS merging. [Default: {SVConfig.MIN_SEQ_SIMILARITY}]")

    # === Group 5: Advanced Options ===
    adv_group = parser.add_argument_group("Advanced Options")
    adv_group.add_argument("--threads", "-t", type=int, default=8, metavar="INT", help="Threads for parallel processing [Default: 8]")
    adv_group.add_argument("--enable-density-filter", action="store_true",
                           help="Enable high-density noise filtering (removes dense clusters of candidates)")
    adv_group.add_argument("--density-window", type=int, default=SVConfig.DENSITY_WINDOW, metavar="INT",
                           help=f"Window size for density filtering (bp). [Default: {SVConfig.DENSITY_WINDOW}]")
    adv_group.add_argument("--max-sv-in-window", type=int, default=SVConfig.MAX_SV_IN_WINDOW, metavar="INT",
                           help=f"Max SVs allowed per density window before suppression. [Default: {SVConfig.MAX_SV_IN_WINDOW}]")
    adv_group.add_argument("--no-retain-sorted", dest="retain_sorted", action="store_false",
                           help="Delete the intermediate queryname-sorted BAM after processing.\n"
                                "By default the .qsort.bam is kept on disk so subsequent runs can\n"
                                "skip the sort step via mtime + integrity reuse check.\n"
                                "Has no effect when --sorted-bam1/2 are provided (user files are never deleted).")
    parser.set_defaults(retain_sorted=True)
    adv_group.add_argument("--force-resort", dest="force_resort", action="store_true",
                           help="Force re-sorting of BAM even if a valid .qsort.bam already exists.\n"
                                "Useful when the cached file is suspected to be stale or from a\n"
                                "different sort run. Has no effect when --sorted-bam1/2 are provided.")
    parser.set_defaults(force_resort=False)
    adv_group.add_argument("--debug", action="store_true", help="Enable verbose DEBUG logging")

    return parser.parse_args(argv)

def update_sv_config(args: argparse.Namespace) -> None:
    """Updates global SVConfig from parsed arguments and calculates derived thresholds."""
    logger.info("Updating global configuration from arguments...")
    
    # 1. Base filter parameters
    if hasattr(args, 'min_sv'): SVConfig.MIN_SV_SIZE = args.min_sv
    if hasattr(args, 'max_sv'): SVConfig.MAX_SV_SIZE = args.max_sv
    if hasattr(args, 'min_mapq'): SVConfig.MIN_MAPQ = args.min_mapq
    if hasattr(args, 'insdel_tolerance'): SVConfig.INSDEL_TOLERANCE = args.insdel_tolerance
    if hasattr(args, 'min_tra_size'): SVConfig.MIN_TRA_SIZE = args.min_tra_size
    if hasattr(args, 'min_csv'): SVConfig.MIN_CSV_SIZE = args.min_csv
    if hasattr(args, 'max_csv'): SVConfig.MAX_CSV_SIZE = args.max_csv

    # 2. Clustering parameters
    if hasattr(args, 'merge_indels_tol'): SVConfig.MERGE_INDEL_TOL = args.merge_indels_tol
    if hasattr(args, 'cluster_pos_tol'): SVConfig.CLUSTER_POS_TOL = args.cluster_pos_tol
    if hasattr(args, 'min_seq_sim'): SVConfig.MIN_SEQ_SIMILARITY = args.min_seq_sim
    
    # 3. Derived length tolerances calculation
    if hasattr(args, 'cluster_len_tol'):
        SVConfig.CLUSTER_LEN_TOL = args.cluster_len_tol
        SVConfig.CLUSTER_del_LEN_TOL = SVConfig._BASE_DEL_LEN_TOL + SVConfig.CLUSTER_LEN_TOL
        SVConfig.CLUSTER_ins_LEN_TOL = SVConfig._BASE_INS_LEN_TOL + SVConfig.CLUSTER_LEN_TOL
    
    # 4. Output controls and flags
    if hasattr(args, 'enable_large_sequence'): SVConfig.ENABLE_LARGE_SEQUENCE = args.enable_large_sequence
    if hasattr(args, 'split_segments'): SVConfig.ENABLE_SPLIT_SEGMENTS = args.split_segments
    if hasattr(args, 'enable_density_filter'): SVConfig.ENABLE_DENSITY_FILTER = args.enable_density_filter
    if hasattr(args, 'density_window'): SVConfig.DENSITY_WINDOW = args.density_window
    if hasattr(args, 'max_sv_in_window'): SVConfig.MAX_SV_IN_WINDOW = args.max_sv_in_window
    
    # 5. Logical integrity checks
    if SVConfig.MIN_SV_SIZE >= SVConfig.MAX_SV_SIZE and SVConfig.MAX_SV_SIZE != -1:
        logger.warning(f"Configuration Warning: MIN_SV ({SVConfig.MIN_SV_SIZE}) is >= MAX_SV ({SVConfig.MAX_SV_SIZE})")

# ==============================================================================
# Sequence & Alignment Utilities
# ==============================================================================

def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a nucleotide sequence."""
    if not seq:
        return ""
    return seq.translate(COMPLEMENT_MAP)[::-1]

def get_target_ref_seq(fasta: pysam.FastaFile, chrom: str, start: int, end: int, strand: str) -> str:
    """
    Safely retrieves the reference sequence for a specified region.
    Returns the reverse complement if the strand is negative.
    """
    try:
        seq = fasta.fetch(chrom, max(0, start), end).upper()
        if not seq: 
            return ""
        return reverse_complement(seq) if strand == '-' else seq
    except (ValueError, KeyError, IndexError):
        return ""

def validate_file(path: str) -> None:
    """Validates if a file exists and is accessible."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file does not exist: {path}")
    if not os.path.isfile(path):
        raise IsADirectoryError(f"Path is a directory, not a file: {path}")

def local_align(ref_seq: str, query_seq: str) -> Tuple[int, int, int, int, float]:
    """
    Performs a local sequence alignment using edlib (NW task=path).
    Returns (ref_start, ref_end, query_start, query_end, normalized_similarity).
    """
    if not ref_seq or not query_seq:
        return 0, 0, 0, 0, 0.0

    result = edlib.align(query_seq, ref_seq, mode="NW", task="path")

    if result.get("editDistance", -1) == -1:
        return 0, 0, 0, 0, 0.0

    edit_dist = result["editDistance"]
    query_start = result.get("queryStart", 0)
    query_end = result.get("queryEnd", 0) + 1 
    ref_start = result.get("targetStart", 0)
    ref_end = result.get("targetEnd", 0) + 1

    align_len = max(ref_end - ref_start, query_end - query_start, 1)
    normalized_sim = max(0.0, 1.0 - (edit_dist / align_len))

    return ref_start, ref_end, query_start, query_end, normalized_sim

# ==============================================================================
# Filtering Algorithms
# ==============================================================================

def filter_candidates_by_density(candidates: List[SVCandidate], window_size: int = 1000, max_count: int = 3) -> List[SVCandidate]:
    """
    Filters high-density SV regions using a sliding window approach (Non-Maximum Suppression variant).
    Retains only the largest SV within a dense window exceeding the max_count limit.
    """
    if not candidates:
        return []

    candidates_by_chr = defaultdict(list)
    for cand in candidates:
        candidates_by_chr[cand.chrom].append(cand)

    final_candidates = []

    for chrom, chrom_cands in candidates_by_chr.items():
        if not chrom_cands:
            continue

        chrom_cands.sort(key=lambda x: x.pos)
        drop_indices: set = set()
        # deque of (sorted_index, candidate) pairs
        window: deque = deque()

        for i, curr_cand in enumerate(chrom_cands):
            while window and (curr_cand.pos - window[0][1].pos > window_size):
                window.popleft()

            window.append((i, curr_cand))

            if len(window) > max_count:
                active = [(idx, c) for idx, c in window if idx not in drop_indices]

                if len(active) > max_count:
                    best_idx = max(active, key=lambda t: t[1].sv_len)[0]
                    for idx, _ in active:
                        if idx != best_idx:
                            drop_indices.add(idx)

        final_candidates.extend(c for i, c in enumerate(chrom_cands) if i not in drop_indices)

    return final_candidates

# ==============================================================================
# I/O & Output Processing
# ==============================================================================

def prepare_output_path(path_str: str, default_filename: str = "results.vcf") -> str:
    """Resolves output directory paths, expands user ~ and ensures parent dirs exist."""
    if not path_str:
        return os.path.abspath(default_filename)

    expanded_path = os.path.abspath(os.path.expanduser(path_str))
    final_path = expanded_path
    
    if os.path.isdir(expanded_path) or path_str.endswith(os.sep):
        final_path = os.path.join(expanded_path, default_filename)
        
    parent_dir = os.path.dirname(final_path)
    if parent_dir and not os.path.exists(parent_dir):
        try:
            os.makedirs(parent_dir, exist_ok=True)
            logger.info(f"Created output directory: {parent_dir}")
        except OSError as e:
            raise OSError(f"Failed to create output directory {parent_dir}: {e}")

    return final_path

def write_vcf(sv_list: List[SVCandidate], output_path: str, header_contigs: Dict[str, int], fasta: pysam.FastaFile) -> None:
    """
    Generates a standard VCFv4.2 file, optimizing memory and disk I/O operations.
    Sorts SVs in biological order (chr1-22, X, Y, M).
    """
    def get_chrom_rank(chrom: str) -> Tuple[int, Any]:
        """Provides sorting weights for chromosome contigs."""
        c = str(chrom).upper()
        if c.startswith('CHR'): c = c[3:]
            
        if c.isdigit(): return (0, int(c))      
        elif c == 'X': return (1, 23)          
        elif c == 'Y': return (1, 24)          
        elif c in ('M', 'MT'): return (1, 25)  
        else: return (2, chrom)                
            
    # Pre-compute chromosome sorting weights for fast O(1) lookups
    chrom_ranks = {chrom: get_chrom_rank(chrom) for chrom in header_contigs.keys()}
    missing_chroms = {sv.chrom for sv in sv_list if sv.chrom not in chrom_ranks}
    for c in missing_chroms:
        chrom_ranks[c] = get_chrom_rank(c)
        
    # Sort structurally using purely deterministic keys
    sv_list.sort(key=lambda x: (chrom_ranks[x.chrom], x.pos))
    
    # Utilize a large file buffer (10MB) to minimize system I/O calls
    try:
        with open(output_path, 'w', buffering=10*1024*1024) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
            f.write("##source=AsmSV\n")
            
            sorted_contigs = sorted(header_contigs.items(), key=lambda item: chrom_ranks[item[0]])
            
            #Batch write Headers and Definitions using writelines
            contig_lines = [f"##contig=<ID={name},length={length}>\n" for name, length in sorted_contigs]
            f.writelines(contig_lines)
                
            alt_defs = [
                ("INS", "Insertion"),
                ("DEL", "Deletion"),
                ("DUP", "Duplication"),
                ("INV", "Inversion"),
                ("BND", "Breakend"),
                ("TRA", "Translocation"),
                ("CSV", "Complex sequence variant"),
                ("COMPLEX", "Complex structural variant"),
            ]
            alt_lines = [f"##ALT=<ID={t},Description=\"{desc}\">\n" for t, desc in alt_defs]
            f.writelines(alt_lines)

            info_defs = [
                # Standard SV fields (VCF 4.2)
                ("SVTYPE", 1, "String", "Type of structural variant (DEL, INS, DUP, INV, BND, TRA, CSV)"),
                ("SVLEN", 1, "Integer", "Signed length difference between REF and ALT alleles; negative for deletions"),
                ("END", 1, "Integer", "End position of the variant on the reference (1-based, inclusive); equals POS for breakends"),
                ("MATEID", 1, "String", "ID of the partner breakend record in a BND pair (VCF 4.2)"),

                # Haplotype-resolved assembly support
                ("SUPPORT_H1", 1, "Integer", "Number of supporting assembly read groups from haplotype 1 after cross-haplotype clustering"),
                ("SUPPORT_H2", 1, "Integer", "Number of supporting assembly read groups from haplotype 2 after cross-haplotype clustering"),
                ("RNAMES", 1, "String", "Comma-separated names of assembly reads (query sequences) supporting this variant"),
                ("SOURCE", 1, "String", "Detection layer: backbone (collinear split-read gap), CIGAR (single-segment CIGAR INDEL), CIGAR_MERGED (merged adjacent CIGAR INDELs), COMPLEX (branch-cluster within backbone gap), COMPLEX_EDGE (edge branch cluster outside backbone)"),

                # Breakend (BND/TRA) mate coordinates
                ("TRA_TYPE", 1, "String", "Breakend scope: INTER_CHR (breakpoints on different chromosomes) or INTRA_CHR (both breakpoints on the same chromosome)"),
                ("CHR2", 1, "String", "Chromosome of the mate breakend"),
                ("POS2", 1, "Integer", "1-based reference position of the mate breakend"),
                ("STRANDS", 1, "String", "Strand orientations at POS and POS2 (first and second character, each + or -), used to construct the BND ALT allele per VCF 4.2"),

                # Complex variant (CSV) decomposition
                ("EVENTTYPE", 1, "String", "Mechanism label on COMPLEX summary records: single branch-cluster subtype, or COMPLEX_<subtype1>_<subtype2>_... for multi-cluster events"),
                ("SUBTYPE", 1, "String", "Mechanistic subclass. Linear: DUP:TANDEM, MULTI_DUP:TANDEM, INS_DEL. Inversion-related: INV, INV-INS, INV-DEL, INV-INS-DEL, InvDUP, InvDUP:TANDEM, InvDUP:PROXIMAL, InvDUP:DISP:INTRA. Dispersed duplication: DUP:PROXIMAL, DUP:DISP:INTRA, DUP:DISP:INTER. Other: LOCAL_RETAIN, COMPLEX"),
                ("DEL_LEN", 1, "Integer", "Length in base pairs of reference sequence deleted or replaced at this locus"),
                ("INS_LEN", 1, "Integer", "Length in base pairs of novel sequence inserted at this locus (query-derived)"),
                ("SEQ", 1, "String", "Inserted sequence; for INS records placed in the ALT column when below length threshold, for CSV(INS_DEL) records placed in INFO"),
                ("BP_DETAILS", 1, "String", "Inversion flanking gap metrics at outer backbone breakpoints: LC/RC = left/right query gaps (bp) between backbone and inverted branch; LF/RF = left/right reference gaps (bp) between outer breakpoints and inverted segment boundaries"),
                ("BP_DEL", 1, "Integer", "Net reference sequence loss (bp) at inversion flanking breakpoints; present only when flanking deletions exceed MIN_SV_SIZE"),
                ("BP_INS", 1, "Integer", "Net inserted sequence length (bp) at inversion flanking breakpoints; present only when flanking insertions exceed MIN_SV_SIZE"),

                # Event linkage and nested structure
                ("EVENT", ".", "String", "Shared event identifier linking related records from one complex rearrangement on an assembly read (BND pairs, span records, summary); comma-separated if merged across clusters"),
                ("ORPHAN_CLUSTERS", 1, "Integer", "Number of branch (orphan) alignment clusters in a COMPLEX summary record"),
                ("MERGED_BLOCKS", 1, "Integer", "Number of adjacent alignment blocks merged into one nested linear SV within a branch network"),

                # Source coordinates and topology trace
                ("SRCREG", 1, "String", "Reference interval of the inserted or duplicated source sequence (format: chr:start-end)"),
                ("SUBREGS", 1, "String", "Comma-separated reference sub-intervals (start-end) when the source region is fragmented across multiple branch alignment blocks; chromosome is that of SRCREG"),
                ("SUB_SVS", 1, "String", "Pipe-separated topology trace on COMPLEX summary records. First entry: backbone,chr,ref_start,ref_end,query_start,query_end,strand; subsequent entries: subtype,chr,ref_start,ref_end,query_start,query_end,strand for each branch cluster"),

                # Inversion and tandem duplication detail
                ("INV_SPAN", 1, "Integer", "Length in base pairs of the inner inverted segment on the reference (inv_end minus inv_start)"),
                ("TANDEM_SOURCE", 1, "String", "Reference interval of the template region for a tandem duplication or short tandem insertion (format: chr:start-end)"),
            ]
            
            
            # 写入 INFO 头部行
            info_lines = [f"##INFO=<ID={i_id},Number={i_num},Type={i_type},Description=\"{i_desc}\">\n" 
                          for i_id, i_num, i_type, i_desc in info_defs]
            f.writelines(info_lines)
            
            filter_defs = [
                ("PASS", "All filters passed"),
                ("COMPLEX", "Multi-cluster (>=2 branch clusters) complex SV summary record; component BND and span records may have FILTER=PASS"),
                ("LOCAL_BND", "Intra-chromosomal breakend with partner distance below MIN_TRA_SIZE; local junction within a complex cluster, not a distal translocation"),
                ("OVERSIZE", "Variant size exceeds the user-defined --max-sv threshold"),
                ("EDGE_COMPLEX", "Complex rearrangement involving branch clusters outside the collinear backbone (read-edge orphans); applies to summary, span, and BND records from that event"),
            ]
            for fid, fdesc in filter_defs:
                f.write(f"##FILTER=<ID={fid},Description=\"{fdesc}\">\n")
            # f.write("##FILTER=<ID=LOCAL_RETAIN,Description=\"Orphan block falls entirely inside the backbone gap (same strand); retained in-gap sequence, not a true duplication\">\n")
            # 写入 FORMAT / 样本列头部
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype: 1|0 (Hap1 only), 0|1 (Hap2 only), 1|1 (both haplotypes), ./. (unknown)\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            
            # Generator-based processing to avoid loading all strings into memory
            type_counter = defaultdict(int)
            type_counter_all=defaultdict(int)
            count = 0
            
            def record_generator():
                nonlocal count
                for sv in sv_list:
                    try:
                        if not sv.id:
                            type_counter[sv.type] += 1
                            current_type_idx = type_counter[sv.type]
                        else:
                            
                            current_type_idx = 0  # unused; sv.id already carries the VCF ID

                        record_tuple = sv.to_vcf_record(
                            current_type_idx,
                            alt_allele=True,
                            reference=fasta
                        )
                        type_counter_all[sv.type] += 1
                        count += 1
                        yield "\t".join(map(str, record_tuple)) + "\n"
                    except Exception as e:
                        logger.error(f"Writing VCF record failed for SV at {sv.chrom}:{sv.pos} | Error: {e}")
                        continue
            
            f.writelines(record_generator())
                
        logger.info(f"Successfully wrote {count} SVs to output file: {output_path}")
        logger.info(f"SV Counts by Type: {dict(type_counter_all)}")

    except IOError as e:
        logger.error(f"Failed to open or write to output VCF file '{output_path}': {e}")