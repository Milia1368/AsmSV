#!/usr/bin/env python3
import os
import sys
import uuid
import time
import shutil
import logging
import subprocess
import pysam
from itertools import groupby
from typing import List, Dict, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED

# Internal module imports
from AsmSV.structures import SVCandidate,BedMaskFilter, filter_candidates_by_mask
from AsmSV.genotype import cluster_and_genotype
from AsmSV.utils import (
    validate_file, filter_candidates_by_density, write_vcf, 
    prepare_output_path, parseArgs, update_sv_config
)
from AsmSV.config import SVConfig
from AsmSV.bam_info import (
    AlignedSegment, split_segments, solve_chain_dp, 
    find_collinear_backbone, resolve_chain_svs
)

logger = logging.getLogger(__name__)

# ==============================================================================
# Core Serialization Module (For safe IPC of BAM records)
# ==============================================================================
def serialize_aligned_segment(read: pysam.AlignedSegment) -> Optional[Dict]:
    """Serializes alignment records while preserving precision (AS & CIGAR)."""
    if read is None or not read.cigartuples:
        return None

    def _get_as_tag(r: pysam.AlignedSegment) -> int:
        try:
            as_score = r.get_tag("AS")
        except KeyError:
            as_score = -1

        # Fallback dynamic scoring mechanism if AS tag is missing or invalid
        if as_score <= 0:
            matched_bases = 0
            small_error_bases = 0  # Sequencing noise (X, and small Indels < 50bp)
            large_sv_events = 0    # True large SV events (>= 50bp)
            
            for op, length in r.cigartuples:
                if op in {0, 7}: 
                    matched_bases += length
                elif op == 8:    
                    small_error_bases += length
                elif op in {1, 2}: 
                    if length < 50:
                        small_error_bases += length  
                    else:
                        large_sv_events += 1         
            
            total_aligned = matched_bases + small_error_bases
            bg_noise_rate = small_error_bases / total_aligned if total_aligned > 0 else 1.0
            
            # Differentiate penalty cost based on mapping quality
            large_cost = 50 if r.mapping_quality >= 50 else 1000
            custom_score = matched_bases * 1.0 - small_error_bases * 2.0 - large_sv_events * large_cost
            
            if (custom_score > 0 and matched_bases >= 500 and 
               (bg_noise_rate < 0.20 or large_sv_events <= 3)):
                as_score = int(custom_score)
            else:
                as_score = -1

        return as_score
    
    final_as_score = _get_as_tag(read)
    if final_as_score < 0:
        return None
        
    return {
        "q_name": read.query_name,
        "ref_name": read.reference_name,
        "ref_start": read.reference_start,
        "ref_end": read.reference_end,
        "cigar": read.cigartuples,
        "mapq": read.mapping_quality,
        "is_reverse": read.is_reverse,
        "flag": read.flag,
        "as_score": final_as_score
    }

def deserialize_aligned_segment(serialized: Dict) -> Optional[AlignedSegment]:
    """Restores the aligned segment from a serialized dictionary."""
    if not serialized: 
        return None
    seg = AlignedSegment(from_dict=serialized.copy())
    seg.is_valid = True
    return seg

def process_single_read_group(read_name: str, serialized_reads_list: List[Dict], query_seq: str, query_len: int) -> List[SVCandidate]:
    """Processes a single read group to extract SV candidates."""
    try:
        raw_segs = []
        for serialized in serialized_reads_list:
            seg = deserialize_aligned_segment(serialized)
            if seg and seg.is_valid:
                raw_segs.append(seg)

        if not raw_segs:
            return []
            
        candidates = []
        
        # (optional)  splitting
        dp_segs = split_segments(raw_segs) if SVConfig.ENABLE_SPLIT_SEGMENTS else raw_segs   
            
        # 1. DP pathfinding and collinear backbone identification
        best_chain = solve_chain_dp(dp_segs, query_len)
        main_chain_idx, main_strand = find_collinear_backbone(best_chain)
        
        main_chain_idx = set(main_chain_idx) if main_chain_idx is not None else set()
            
        # 2. Resolve structural variations from the alignment chain
        splits_svs = resolve_chain_svs(best_chain, main_strand, main_chain_idx, read_name, query_seq)
        candidates.extend(splits_svs)
        
        # Free memory
        del raw_segs
        
        if not best_chain:
            return []

        # 3. parsing intra-alignment linear variations (CIGAR SVs)
        for seg in best_chain:
            raw_svs = seg.parse_intra_alignment_sv(
                query_seq=query_seq, 
                min_size=SVConfig.MIN_SV_SIZE, 
                merge_gap_tol=SVConfig.MERGE_INDEL_TOL
            )
            for s in raw_svs:
                candidates.append(SVCandidate(s['type'], s['chr'], s['pos'], s['len'], s['read'], s['info']))

        return candidates
    except Exception as e:
        logger.error(f"Error processing read {read_name}: {e}", exc_info=logger.isEnabledFor(logging.DEBUG))
        return [] 
    

def worker_wrapper(chunk_of_read_groups: List[Tuple]) -> List[SVCandidate]:
    """Wrapper function to process chunks of read groups in worker processes."""
    results = []
    for name, serialized_reads, query_seq, query_len in chunk_of_read_groups:
        res = process_single_read_group(name, serialized_reads, query_seq, query_len)
        results.extend(res)
    return results

# ==============================================================================
# BAM Data Pump (Read extraction & preprocessing)
# ==============================================================================
def bam_generator(bam_path: str):
    """Yields serialized read groups from a name-sorted BAM file."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for qname, reads_iter in groupby(bam, key=lambda x: x.query_name):
            filtered_reads = [r for r in reads_iter if r.mapping_quality >= SVConfig.MIN_MAPQ]
            if not filtered_reads:
                continue
                
            query_seq, query_len = None, None
            # Filter out secondary and supplementary alignments
            primary_alignments = [aln for aln in filtered_reads if not (aln.flag & 0x900)] 
            
            # Enforce get_forward_sequence() to ensure uniform 5'->3' sequences
            if primary_alignments:
                best_primary = max(primary_alignments, key=lambda x: x.mapping_quality)
                query_seq = best_primary.get_forward_sequence()
                query_len = best_primary.infer_query_length()
            else:
                for aln in filtered_reads:
                    fwd_seq = aln.get_forward_sequence()
                    if fwd_seq is not None:
                        query_seq = fwd_seq
                        query_len = aln.infer_query_length()
                        break

            # Discard read if forward sequence cannot be obtained
            if not query_seq or not query_len:
                continue

            serialized_reads = [s for r in filtered_reads if (s := serialize_aligned_segment(r)) is not None]
            if serialized_reads:
                yield qname, serialized_reads, query_seq, query_len

def sort_bam_by_name(input_bam: str, threads: int = 4, reuse_sorted: bool = True) -> Tuple[str, bool]:
    """Sorts the BAM file by query name."""
    sorted_bam = input_bam.replace(".bam", ".qsort.bam")

    if reuse_sorted and os.path.exists(sorted_bam):
        try:
            if os.path.getmtime(sorted_bam) >= os.path.getmtime(input_bam):
                with pysam.AlignmentFile(sorted_bam, "rb") as bam:
                    _ = bam.header  
                logger.info(f"Reusing qsort BAM: {sorted_bam}")
                return sorted_bam, False
        except Exception:
            pass

    logger.info(f"Checking integrity of input BAM: {input_bam}")

    if shutil.which("samtools"):
        check_result = subprocess.run(
            ["samtools", "quickcheck", "-v", input_bam],
            capture_output=True, text=True,
        )
        if check_result.returncode != 0:
            logger.error(f"INPUT FILE CORRUPTED: The file {input_bam} is truncated or damaged!")
            logger.error(f"Reason: {check_result.stdout.strip() or check_result.stderr.strip()}")
            raise ValueError(f"Corrupted input BAM file: {input_bam}")
    else:
        logger.warning("samtools not found in PATH, skipping BAM integrity check.")

    file_size_gb = os.path.getsize(input_bam) / (1024 ** 3)
    # if file_size_gb < 1.0:
    #     logger.warning(f"WARNING: The input BAM file is unusually small ({file_size_gb:.2f} GB).")

    if os.path.exists(sorted_bam):
        try: os.remove(sorted_bam)
        except OSError: pass
            
    # Clean up invalid indices
    for ext in (".bai", ".csi"):
        idx_path = sorted_bam + ext
        if os.path.exists(idx_path):
            try: os.remove(idx_path)
            except OSError: pass

    temp_prefix = sorted_bam + ".tmp." + str(uuid.uuid4().hex[:8])
    logger.info(f"Sorting BAM by query name: {input_bam} -> {sorted_bam}")

    try:
        pysam.sort("-n", "-@", str(threads), "-o", sorted_bam, "-T", temp_prefix, input_bam)
    except Exception as e:
        logger.error(f"Sorting failed: {e}")
        if os.path.exists(sorted_bam):
            os.remove(sorted_bam)
        raise

    return sorted_bam, True

# ==============================================================================
# Core Parallel Execution Pipeline
# ==============================================================================
def main_pipeline(bam_path: str,label:str, mask_filter: BedMaskFilter,threads: int = 8) -> Tuple[List[SVCandidate], Dict]:
    """Manages multiprocessing execution for processing SVs from a BAM file."""
    header_dict = {}
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            bam_header_sq = bam.header.to_dict().get('SQ', [])
            header_dict = {d['SN']: d['LN'] for d in bam_header_sq}
    except Exception as e:
        logger.error(f"Failed to read BAM header from {bam_path}: {e}")
        return [], {}

    qsort_bam, created_now = sort_bam_by_name(bam_path, threads=threads, reuse_sorted=True)
    
    # Configure process pool and control IPC throughput
    executor = ProcessPoolExecutor(max_workers=threads)
    chunk_size = 15
    current_chunk = []
    pending_futures = set()
    all_svs = []
    
    # Restrict queued tasks in memory to prevent OOM errors
    max_in_flight = max(1, threads * 3) 
    logger.info(f"Start processing reads from {qsort_bam}...")

    try:
        for payload in bam_generator(qsort_bam):
            current_chunk.append(payload)
            
            if len(current_chunk) >= chunk_size:
                fut = executor.submit(worker_wrapper, current_chunk)
                pending_futures.add(fut)
                current_chunk = []

                if len(pending_futures) >= max_in_flight:
                    done_futures, pending_futures = wait(pending_futures, return_when=FIRST_COMPLETED)
                    for f in done_futures:
                        try:
                            all_svs.extend(f.result())
                        except Exception as e:
                            logger.error(f"Worker future failed: {e}")

        if current_chunk:
            pending_futures.add(executor.submit(worker_wrapper, current_chunk))

        for future in as_completed(pending_futures):
            try:
                all_svs.extend(future.result())
            except Exception as e:
                logger.error(f"Worker future failed: {e}")

    finally:
        executor.shutdown(wait=True)
        # Only delete qsort.bam if it was generated during this execution
        if created_now and qsort_bam != bam_path and os.path.exists(qsort_bam):
            try: os.remove(qsort_bam)
            except OSError: pass

    # 区域屏蔽过滤
    if mask_filter and mask_filter.is_active:
        all_candidates = filter_candidates_by_mask(all_candidates, mask_filter)
        logger.info(f"{label} After Masking: {len(all_candidates)}")

    logger.info(f"Total SV candidates found in {os.path.basename(bam_path)}: {len(all_svs)}")

    if SVConfig.ENABLE_DENSITY_FILTER :
        filtered_svs = filter_candidates_by_density(
            all_svs, 
            window_size=SVConfig.DENSITY_WINDOW, 
            max_count=SVConfig.MAX_SV_IN_WINDOW
        )
        return filtered_svs, header_dict
    
    return all_svs, header_dict

# ==============================================================================
# Main Controller & Entry Point
# ==============================================================================
def main_ctrl(args, argv):
    """Orchestrates configuration, validation, and pipeline execution."""
    update_sv_config(args)
    cmd_line = "asmSV " + (" ".join(argv) if isinstance(argv, list) else str(argv))

    logger.info("[Command Line]\n" f"{cmd_line}")
    logger.info("================ Execution Configuration ================")
    logger.info(f" - Threads         : {args.threads}")
    logger.info(f" - Min/Max SV Size : {SVConfig.MIN_SV_SIZE} / {SVConfig.MAX_SV_SIZE}")
    logger.info(f" - Min MapQ        : {SVConfig.MIN_MAPQ}")
    logger.info(f" - Min TRA Size    : {SVConfig.MIN_TRA_SIZE}")
    logger.info(f" - Density Filter  : {'Enabled' if SVConfig.ENABLE_DENSITY_FILTER else 'Disabled'}")
    logger.info(f" - Input Ref       : {args.ref}")
    logger.info(f" - Input Hap1 BAM  : {args.bam1}")
    logger.info(f" - Input Hap2 BAM  : {args.bam2}")
    if getattr(args, 'mask', None):
        logger.info(f" - Mask BED        : {args.mask}")
    logger.info("=========================================================")

    try:
        final_output_path = prepare_output_path(args.output, "results.vcf")
        logger.info(f"Output VCF will be saved to: {final_output_path}")
    except Exception as e:
        logger.error(f"Output path error: {e}")
        sys.exit(1)

    try:
        validate_file(args.bam1)
        validate_file(args.bam2)
        validate_file(args.ref)
    except Exception as e:
        logger.error(f"Input validation failed: {e}")
        sys.exit(1)

    ref_fasta = pysam.FastaFile(args.ref)
   
    mask_filter = BedMaskFilter(getattr(args, 'mask', None))
    start_time = time.time()
    
    try:
        # Process haplotypes concurrently
        svs_h1, header_dict = main_pipeline(args.bam1, "hap1",mask_filter,args.threads)
        svs_h2, _ = main_pipeline(args.bam2, "hap2",mask_filter,args.threads)

        logger.info(">>> Global Clustering and Genotyping...")
        final_svs = cluster_and_genotype(svs_h1, svs_h2)

        write_vcf(final_svs, final_output_path, header_dict, ref_fasta)

        elapsed = time.time() - start_time
        logger.info(f"Pipeline finished successfully in {elapsed:.2f}s")

    except Exception as e:
        logger.error(f"Critical failure during execution: {e}", exc_info=True)
        sys.exit(1)

def setupLogging(debug=False, log_file="asmSV.log"):
    """Configures system logging outputs."""
    log_level = logging.DEBUG if debug else logging.INFO
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    formatter = logging.Formatter(log_format)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # File handler
    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

def run(argv):
    """Application entry function."""
    args = parseArgs(argv)
    setupLogging(args.debug)
    starttime = time.time()
    main_ctrl(args, argv)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))

def main():
    """
    Parameterless entry point for setuptools console_scripts.
    It automatically reads command-line arguments from sys.argv.
    """
    sys.exit(run(sys.argv[1:]))

if __name__ == '__main__':
    main()