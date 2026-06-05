import logging
import edlib
from typing import List, Callable, Dict, Set, Optional
from collections import defaultdict

from AsmSV.structures import SVCandidate
from AsmSV.config import SVConfig

logger = logging.getLogger(__name__)

# ==============================================================================
# Base Utilities
# ==============================================================================

def calculate_seq_similarity(seq1: str, seq2: str) -> float:
    """
    Calculates sequence similarity tailored for assembly data.
    Retains physical edge trimming to mitigate aligner wobble but tightens internal edit tolerance.
    """
    if not seq1 or not seq2:
        return 0.0

    len1, len2 = len(seq1), len(seq2)

    # Mechanism 1: Flank Trimming (mitigates wobble in tandem repeats)
    trim_len = 3 
    if len1 > 3 * trim_len and len2 > 3 * trim_len:
        seq1_core = seq1[trim_len:-trim_len]
        seq2_core = seq2[trim_len:-trim_len]
    else:
        seq1_core, seq2_core = seq1, seq2

    core_max_len = max(len(seq1_core), len(seq2_core))
    align_res = edlib.align(seq1_core, seq2_core, mode="NW", task="distance")
    
    if align_res["editDistance"] == -1:
        return 0.0
        
    ed = align_res["editDistance"]

    # Mechanism 2: Base Edit Tolerance (allows minimal polishing residual errors)
    BASE_EDIT_TOLERANCE = 0
    if ed <= BASE_EDIT_TOLERANCE:
        return 1.0
        
    adjusted_ed = ed - BASE_EDIT_TOLERANCE
    return 1.0 - (adjusted_ed / core_max_len)

# ==============================================================================
# SV Comparators
# ==============================================================================

def _match_span_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic for DEL, DUP, and INV.
    Prioritizes positional and length proximity over sequence alignment for speed.
    """
    if s1.chrom != s2.chrom: 
        return False
    
    if abs(s1.pos - s2.pos) > SVConfig.CLUSTER_POS_TOL:
        return False

    len_diff = abs(s1.sv_len - s2.sv_len)
    if len_diff < SVConfig.CLUSTER_del_LEN_TOL:
        return True
        
    return False

def _match_ins_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic for INS.
    Requires strict positional proximity and sequence similarity.
    """
    if s1.chrom != s2.chrom: 
        return False
    
    if abs(s1.pos - s2.pos) > SVConfig.CLUSTER_POS_TOL:
        return False

    if abs(s1.sv_len - s2.sv_len) > SVConfig.CLUSTER_ins_LEN_TOL:
        return False

    seq1 = s1.info.get("SEQ", "")
    seq2 = s2.info.get("SEQ", "")

    if seq1 and seq2:
        sim = calculate_seq_similarity(seq1, seq2)
        return sim >= SVConfig.MIN_SEQ_SIMILARITY
    
    # Fallback if sequence is missing
    return True

def _match_bnd_tra_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic for BND and TRA (Translocations).
    Requires both breakpoints to be proximally matched.
    """
    if s1.chrom != s2.chrom: 
        return False
    
    if abs(s1.pos - s2.pos) > SVConfig.CLUSTER_POS_TOL: 
        return False
        
    chr2_1 = s1.info.get("CHR2")
    chr2_2 = s2.info.get("CHR2")
    if chr2_1 != chr2_2: 
        return False
        
    pos2_1 = s1.info.get("POS2", -1) 
    pos2_2 = s2.info.get("POS2", -1)
    if abs(pos2_1 - pos2_2) > SVConfig.CLUSTER_POS_TOL: 
        return False
        
    return True
def _match_csv_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic specifically for CSV (Complex Sequence Replacement: INS + DEL).
    Checks positional boundaries (POS and END) and compares INS_LEN and DEL_LEN components separately.
    """
    # 1. 染色体必须相同
    if s1.chrom != s2.chrom: 
        return False
    
    # 2. 起点位置 (POS) 必须相近
    if abs(s1.pos - s2.pos) > SVConfig.CLUSTER_POS_TOL:
        return False
        
    # 3. 终点位置 (END) 必须相近
    # 这非常关键，对于 CSV 来说，END 锚定了 Reference 上被替换序列的物理右边界
    if abs(s1.end - s2.end) > SVConfig.CLUSTER_POS_TOL:
        return False

    # 4. 提取独立的插入和缺失长度
    ins_len1 = int(s1.info.get("INS_LEN", 0))
    del_len1 = int(s1.info.get("DEL_LEN", 0))
    
    ins_len2 = int(s2.info.get("INS_LEN", 0))
    del_len2 = int(s2.info.get("DEL_LEN", 0))

    diff_ins = abs(ins_len1 - ins_len2)
    diff_del = abs(del_len1 - del_len2)
    
    # 获取容差设定 (这里可以复用 del 的容差，如果没专门配 ins 的容差，就 fallback 回 del_len_tol)
    ins_tol = SVConfig.CLUSTER_ins_LEN_TOL
    del_tol = SVConfig.CLUSTER_del_LEN_TOL

    # 5. 核心阈值判断：插入部分和缺失部分的长度差异必须同时满足要求
    # 额外加入动态比例判断：如果绝对差异略微超出阈值，但差异占比不到该维度的 20%，也可以认为是同一个复杂变异的测序摆动
    ins_match = (diff_ins <= ins_tol) or (diff_ins <= max(ins_len1, ins_len2) * 0.1)
    del_match = (diff_del <= del_tol) or (diff_del <= max(del_len1, del_len2) * 0.1)

    if ins_match and del_match:
        return True
        
    return False
# ==============================================================================
# Core Clustering Engine
# ==============================================================================

def _core_cluster_logic(sv_list: List[SVCandidate], match_func: Callable) -> List[SVCandidate]:
    """
    Generic linear clustering logic.
    Assumes sv_list is already sorted by (chrom, pos).
    """
    if not sv_list:
        return []
        
    results = []
    current_cluster = [sv_list[0]]
    
    for curr_sv in sv_list[1:]:
        rep_sv = current_cluster[0] # Compare with cluster representative
        
        if match_func(rep_sv, curr_sv):
            current_cluster.append(curr_sv)
        else:
            merged = _merge_candidates(current_cluster)
            if merged: 
                results.append(merged)
            current_cluster = [curr_sv]
            
    if current_cluster:
        merged = _merge_candidates(current_cluster)
        if merged: 
            results.append(merged)
        
    return results

def _merge_candidates(cluster: List[SVCandidate]) -> Optional[SVCandidate]:
    """
    Merges a cluster of SV candidates into a single Genotyped SV.
    """
    if not cluster: 
        return None
    if len(cluster) == 1:
        sv = cluster[0]
        has_h1 = "Hap1" in sv.sources
        has_h2 = "Hap2" in sv.sources
        
        if has_h1 and has_h2: 
            gt = "1|1"
        elif has_h1: 
            gt = "1|0"
        elif has_h2: 
            gt = "0|1"
        else: 
            gt = "./."
            
        sv.genotype = gt
        sv.info["SUPPORT_H1"] = 1 if has_h1 else 0
        sv.info["SUPPORT_H2"] = 1 if has_h2 else 0
        return sv
    cluster.sort(key=lambda x: x.sv_len)
    rep = cluster[len(cluster) // 2]
    
    sources = set()
    all_read_names = set()
    h1_cnt, h2_cnt = 0, 0
    merged_events=set()
    merged_pvents=set()
    for sv in cluster:
        if hasattr(sv, 'read_names'):
            all_read_names.update(sv.read_names)
        if "EVENT" in sv.info:
            merged_events.add(sv.info["EVENT"])
        if "PEVENT" in sv.info:
            merged_pvents.add(sv.info["PEVENT"])
        sources.update(sv.sources)
        if "Hap1" in sv.sources: h1_cnt += 1
        if "Hap2" in sv.sources: h2_cnt += 1
            
    # Resolve Phased Genotype
    has_h1, has_h2 = h1_cnt > 0, h2_cnt > 0
    if has_h1 and has_h2: gt = "1|1"
    elif has_h1: gt = "1|0"
    elif has_h2: gt = "0|1"
    else: gt = "./." 
    
    new_sv = SVCandidate(
        sv_type=rep.raw_type,
        chrom=rep.chrom,
        pos=rep.pos,
        sv_len=rep.sv_len,
        read_name="", # Managed via set update below
        info=rep.info.copy(),
        end=rep.end
    )
    new_sv.sources = sources
    new_sv.read_names = all_read_names
    new_sv.genotype = gt
    new_sv.info["SUPPORT_H1"] = h1_cnt
    new_sv.info["SUPPORT_H2"] = h2_cnt

    # 3-B: 保留 mateid（BND 配对 ID）：若簇内所有记录的 mateid 一致则继承
    mateids = {sv.mateid for sv in cluster if sv.mateid}
    if len(mateids) == 1:
        new_sv.mateid = mateids.pop()
    if len(merged_events)!=0:
        # 加上 sorted 可以保证每次输出的顺序一致，不会乱跳
        new_sv.info["EVENT"] = ",".join(sorted(merged_events))
    if len(merged_pvents)!=0:
        new_sv.info["PEVENT"] = ",".join(sorted(merged_pvents))
    # Retain the longest sequence as the representative REF/ALT for INS
    if rep.type == "INS":
        longest_seq_sv = max(cluster, key=lambda x: len(x.info.get("SEQ", "")))
        if longest_seq_sv.info.get("SEQ"):
            new_sv.info["SEQ"] = longest_seq_sv.info.get("SEQ")
            new_sv.sv_len = longest_seq_sv.sv_len 
            
    return new_sv

# ==============================================================================
# Main Pipeline Entry
# ==============================================================================

def cluster_and_genotype(hap1_svs: List[SVCandidate], hap2_svs: List[SVCandidate]) -> List[SVCandidate]:
    """
    Main entry point: segregates complex events, dispatches by SV type, sorts, and clusters.
    """
    # 1. Preprocessing: Assign haplotype genotypes and sources
    for sv in hap1_svs:
        sv.sources = {"Hap1"}
        sv.genotype = "1|0"
    for sv in hap2_svs:
        sv.sources = {"Hap2"}
        sv.genotype = "0|1"
        
    raw_svs = hap1_svs + hap2_svs
    if not raw_svs: 
        return []

    final_results = []
    clusterable_buckets = defaultdict(list)
    
    # 2. Segregate events: FILTER=COMPLEX summary records pass through directly;
    #    all others enter clustering buckets
    for sv in raw_svs:
        if sv.filter == "COMPLEX":
            # 3-A: summary records from [A3] — not meaningful to cluster
            final_results.append(sv)
        elif sv.type == "CSV":
            subtype = sv.info.get("SUBTYPE", "UNKNOWN")
            clusterable_buckets[f"CSV_{subtype}"].append(sv)
        else:
            clusterable_buckets[sv.type].append(sv)

    common_sort_key = lambda x: (x.chrom, x.pos, x.sv_len)
    
    # A. Process DEL, DUP, INV
    for sv_type in ["DEL", "DUP", "INV"]:
        if candidates := clusterable_buckets[sv_type]:
            candidates.sort(key=common_sort_key)
            final_results.extend(_core_cluster_logic(candidates, _match_span_sv))
        
    # B. Process INS
    if ins_candidates := clusterable_buckets["INS"]:
        ins_candidates.sort(key=common_sort_key)
        final_results.extend(_core_cluster_logic(ins_candidates, _match_ins_sv))
        
    # C. Process BND and TRA
    for sv_type in ["BND", "TRA"]:
        if bnd_tra_candidates := clusterable_buckets[sv_type]:
            bnd_tra_candidates.sort(key=lambda x: (x.chrom, x.pos, x.info.get("CHR2", ""), x.info.get("POS2", 0)))
            final_results.extend(_core_cluster_logic(bnd_tra_candidates, _match_bnd_tra_sv))
    # CSV routing:
    #   INS_DEL / INV-* → _match_csv_sv (checks pos+end+ins_len+del_len separately)
    #   InvDUP:* / DUP:* dispersed/proximal → _match_span_sv (pos+sv_len)
    # Note: DUP:TANDEM creates SVTYPE=INS or DUP (not CSV), so no CSV_DUP:TANDEM bucket.
    # Note: UNPAIRED_INV_BP creates SVTYPE=BND, so no CSV_UNPAIRED_INV_BP bucket.
    CSV_INS_DEL_TYPES = {"INS_DEL", "INV-INS", "INV-DEL", "INV-INS-DEL"}
    CSV_SPAN_TYPES = {
        "InvDUP:TANDEM", "InvDUP:PROXIMAL", "InvDUP:DISP:INTRA", "InvDUP:DISP:INTER",
        "DUP:PROXIMAL", "DUP:DISP:INTRA", "DUP:DISP:INTER",
    }

    for sv_type in sorted(CSV_INS_DEL_TYPES | CSV_SPAN_TYPES):
        if csv_candidates := clusterable_buckets.get(f"CSV_{sv_type}"):
            csv_candidates.sort(key=common_sort_key)
            if sv_type in CSV_INS_DEL_TYPES:
                final_results.extend(_core_cluster_logic(csv_candidates, _match_csv_sv))
            else:
                final_results.extend(_core_cluster_logic(csv_candidates, _match_span_sv))


    logger.info(f"Genotyping Summary: Input={len(raw_svs)}, Output={len(final_results)}")
    final_results.sort(key=lambda x: (x.chrom, x.pos))
    
    return final_results