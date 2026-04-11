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
        description="Advanced Long-read Assembly SV Detector (v1.0)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # === Group 1: Input/Output (Required) ===
    io_group = parser.add_argument_group("Input/Output Arguments")
    io_group.add_argument("--bam1", required=True, metavar="FILE", help="Input BAM for Hap1 (sorted & indexed)")
    io_group.add_argument("--bam2", required=True, metavar="FILE", help="Input BAM for Hap2 (sorted & indexed)")
    io_group.add_argument("--ref", required=True, metavar="FILE", help="Reference Genome FASTA (indexed)")
    io_group.add_argument("--output", "-o", metavar="PATH", default="./results.vcf",
                          help="Output VCF path. If a directory is given, defaults to results.vcf inside it. [Default: ./results.vcf]")
    
    # === Group 2: VCF Format Options ===
    fmt_group = parser.add_argument_group("VCF Format Options")
    fmt_group.add_argument("--max-alleles-len", type=int, default=5000, metavar="INT",
                           help="Maximum sequence length to write in ALT column.\n"
                                "SVs longer than this will use symbolic alleles (e.g. <DEL>).[Default: 5000]")
    # Logic: Symbolic alleles are enabled by default. '--no-symbolic' disables it.
    fmt_group.add_argument("--no-symbolic", dest="symbolic_alleles", action='store_false',
                           help="Disable symbolic alleles (e.g. <DEL>) entirely and force sequence output.\n"
                                "WARNING: This may create huge VCF files.")
    parser.set_defaults(symbolic_alleles=True)

    # === Group 3: Variant Filtering Parameters ===
    filter_group = parser.add_argument_group("Filter Parameters")
    filter_group.add_argument("--min-sv", type=int, default=SVConfig.MIN_SV_SIZE, metavar="INT", 
                              help=f"Minimum SV size [Default: {SVConfig.MIN_SV_SIZE}]")
    filter_group.add_argument("--max-sv", type=int, default=SVConfig.MAX_SV_SIZE, metavar="INT",
                              help=f"Maximum SV size [Default: {SVConfig.MAX_SV_SIZE}]. SVs larger than this value will be skipped.\n"
                                   f"Set to -1 to report all variants without upper bound filtering.")
    filter_group.add_argument("--min-mapq", type=int, default=SVConfig.MIN_MAPQ, metavar="INT",
                              help=f"Minimum mapping quality [Default: {SVConfig.MIN_MAPQ}]")
    filter_group.add_argument("--insdel-tolerance", type=int, default=SVConfig.INSDEL_TOLERANCE, metavar="INT",
                              help="Tolerance threshold (bp) for calling split-read INS/DELs.\n"
                                   f"[Default: {SVConfig.INSDEL_TOLERANCE}]")
    filter_group.add_argument("--min-tra-size", type=int, default=SVConfig.MIN_TRA_SIZE, metavar="INT",
                              help="Threshold size to classify a large gap as a Translocation (TRA/BND).\n"
                                   f"[Default: {SVConfig.MIN_TRA_SIZE}]")
    filter_group.add_argument("--split-segments", action="store_true", 
                              help="Enable segment splitting/refinement before DP chaining. [Default: off]")
    
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
    if hasattr(args, 'max_alleles_len'): SVConfig.MAX_ALLELE_LEN = args.max_alleles_len
    if hasattr(args, 'symbolic_alleles'): SVConfig.ENABLE_ALLELE = args.symbolic_alleles
    if hasattr(args, 'split_segments'): SVConfig.ENABLE_SPLIT_SEGMENTS = args.split_segments
    if hasattr(args, 'enable_density_filter'): SVConfig.ENABLE_DENSITY_FILTER = args.enable_density_filter
    
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

def retrieve_other_alignments(main_alignment: pysam.AlignedSegment, bam: pysam.AlignmentFile) -> List[pysam.AlignedSegment]:
    """Extracts supplementary alignments from the SA tag of a read."""
    try:
        sa_tag = main_alignment.get_tag("SA")
    except KeyError:
        return []

    other_alignments = []
    for element in sa_tag.split(";"):
        if not element:
            continue
            
        fields = element.split(",")
        if len(fields) < 6:
            continue
            
        rname, pos, strand, cigar, mapq, nm = fields[:6]
        
        try:
            a = pysam.AlignedSegment()
            a.query_name = main_alignment.query_name
            a.flag = 2048 | (16 if strand == "-" else 0)
            
            tid = bam.get_tid(rname)
            if tid < 0: continue 
            a.reference_id = tid
            a.reference_start = int(pos) - 1  
            a.mapping_quality = min(int(mapq), 255)
            a.cigarstring = cigar
            
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            a.query_sequence = "" 
            a.set_tags([("NM", int(nm), "i")])
            
            other_alignments.append(a)
        except (ValueError, IndexError, TypeError):
            continue
            
    return other_alignments

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
        drop_ids = set()
        window_cands = deque()
        
        for curr_cand in chrom_cands:
            while window_cands and (curr_cand.pos - window_cands[0].pos > window_size):
                window_cands.popleft()
            
            window_cands.append(curr_cand)
            
            if len(window_cands) > max_count:
                active_cands = [c for c in window_cands if id(c) not in drop_ids]
                
                if len(active_cands) > max_count:
                    best_cand = max(active_cands, key=lambda x: x.sv_len)
                    for c in active_cands:
                        if c is not best_cand:
                            drop_ids.add(id(c))
        
        for cand in chrom_cands:
            if id(cand) not in drop_ids:
                final_candidates.append(cand)

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
            f.write("##source=AssembleSV_Refine_Tool\n")
            
            sorted_contigs = sorted(header_contigs.items(), key=lambda item: chrom_ranks[item[0]])
            
            #Batch write Headers and Definitions using writelines
            contig_lines = [f"##contig=<ID={name},length={length}>\n" for name, length in sorted_contigs]
            f.writelines(contig_lines)
                
            alt_lines = [f"##ALT=<ID={t},Description=\"Structural Variant {t}\">\n" 
                         for t in ["INS", "DEL", "DUP", "INV", "TRA", "BND", "COMPLEX"]]
            f.writelines(alt_lines)
            
            info_defs = [
                ("SVTYPE", 1, "String", "Type of structural variant"),
                ("SVLEN", ".", "Integer", "Difference in length between REF and ALT alleles"),
                ("END", 1, "Integer", "End position of the variant described in this record."),
                ("SUPPORT_H1", 1, "Integer", "Number of reads/haplotypes from hap1 supporting the variant"),
                ("SUPPORT_H2", 1, "Integer", "Number of reads/haplotypes from hap2 supporting the variant"),
                ("RNAMES", ".", "String", "Names of supporting reads"),
                ("SOURCE", 1, "String", "Algorithm source: Split_Read or CIGAR"),
                ("STRANDS", 1, "String", "Strand orientation of the structural variant (e.g., +-, -+)"),
                ("SUBTYPE", 1, "String", "Detailed classification of the complex variant (e.g., INV-INS-DEL, InvTRA, InvDUP)"),
                ("TRA_TYPE", 1, "String", "Translocation type: INTER_CHR (diff chroms) or INTRA_CHR (same chrom dist > limit)"),
                ("CHR2", 1, "String", "Chromosome for BND/TRA connection"),
                ("POS2", 1, "Integer", "Position for BND/TRA connection"),
                ("EVENT", 1, "String", "ID of the primary event associated with this breakend or variant"),
                ("PEVENT", 1, "String", "Parent event ID for nested sub-variants in complex SVs"),
                ("REPLACED_REF_LEN", 1, "Integer", "Length of the reference sequence deleted/replaced by a COMPLEX event"),
                ("ORPHAN_CLUSTERS", 1, "Integer", "Number of unmapped or orphan read clusters forming the COMPLEX variant"),
                ("SOURCE_CHR", 1, "String", "Source chromosome of the duplicated or translocated sequence"),
                ("SOURCE_POS", 1, "Integer", "Source start coordinate of the duplicated or translocated sequence"),
                ("SOURCE_END", 1, "Integer", "Source end coordinate of the duplicated or translocated sequence"),
                ("MERGED_BLOCKS", 1, "Integer", "Number of tandem alignment blocks merged into this single SV event"),
                ("COPIES", 1, "Integer", "Estimated copy number of the tandemly amplified sequence within the SV"),
                ("DEST_CHR", 1, "String", "Destination chromosome where the sequence entity is inserted or attached"),
                ("DEST_POS", 1, "Integer", "Destination coordinate where the sequence entity is inserted or attached"),
                ("INV_SPAN", 1, "Integer", "Length of the inverted segment on the reference genome"),
                ("FLANK_SPAN", 1, "Integer", "Span of the flanking normal sequences; negative value indicates Target Site Duplication (TSD)"),
                ("LEFT_DEL", 1, "Integer", "Length of sequence deletion at the left breakpoint of a complex inversion"),
                ("RIGHT_DEL", 1, "Integer", "Length of sequence deletion at the right breakpoint of a complex inversion"),
                ("LEFT_INS", 1, "Integer", "Length of sequence insertion at the left breakpoint of a complex inversion"),
                ("RIGHT_INS", 1, "Integer", "Length of sequence insertion at the right breakpoint of a complex inversion")
            ]
            info_lines = [f"##INFO=<ID={i_id},Number={i_num},Type={i_type},Description=\"{i_desc}\">\n" 
                          for i_id, i_num, i_type, i_desc in info_defs]
            f.writelines(info_lines)
                
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            
            # Generator-based processing to avoid loading all strings into memory
            type_counter = defaultdict(int)
            count = 0
            
            def record_generator():
                nonlocal count
                for sv in sv_list:  
                    try:
                        type_counter[sv.type] += 1
                        current_type_idx = type_counter[sv.type]
                        
                        record_tuple = sv.to_vcf_record(
                            current_type_idx,  
                            alt_allele=SVConfig.ENABLE_ALLELE, 
                            reference=fasta
                        )
                        count += 1
                        yield "\t".join(map(str, record_tuple)) + "\n"
                    except Exception as e:
                        logger.error(f"Writing VCF record failed for SV at {sv.chrom}:{sv.pos} | Error: {e}")
                        continue
            
            f.writelines(record_generator())
                
        logger.info(f"Successfully wrote {count} SVs to output file: {output_path}")
        logger.info(f"SV Counts by Type: {dict(type_counter)}")

    except IOError as e:
        logger.error(f"Failed to open or write to output VCF file '{output_path}': {e}")