# ==============================================================================
# Configuration & Static Constants
# ==============================================================================
from typing import List, Dict, Tuple, Set, Optional, Any
from dataclasses import dataclass, field as dc_field
from AsmSV.config import SVConfig
from AsmSV.structures import SVCandidate

# --- Base Alignment Parameters ---
MIN_ANCHOR_LEN = 10     
ANCHOR_LOOKAHEAD = 3    
SPLIT_SCORE_RATIO = 1.1 
MATCH_SCORE_DEFAULT = 1 
SCORE_MATCH = 1
SCORE_MISMATCH = -4 

# --- Split and Tolerance Parameters ---
REFINE_TOL = 500
HALF_REFINE_TOL = 200

# --- DP Pathfinding and Topology Parameters ---
PENALTY_INV = 200.0         
PENALTY_DUP = 60.0           
PENALTY_BND = 80.0           
MAX_PENALTY_SV = 500.0       # SV penalty cap
SV_JUMP_PENALTY = 50.0       

LINEAR_GAP_RATE = 0.1        # Penalty coefficient for linear gaps
MAX_LINEAR_GAP = 10000       # Distances exceeding this are considered SVs
UNMAPPED_RATE = 0.5          # Penalty for unaligned regions

OVERLAP_PENALTY_RATE = 1.0   # Overlap offset ratio
MAX_ALLOWED_OVERLAP = 10000  # Max permitted overlap threshold (avoids aligner hallucination)
MIN_AS_PER_BP = 0.4          

# Static DNA reverse complement translation table (Performance optimization)
RC_TRANS = str.maketrans("ATCGNatcgn", "TAGCNtagcn")

# ==============================================================================
# Core Data Structure: AlignedSegment
# ==============================================================================
class AlignedSegment:
    __slots__ = [
        'q_name', 'flag', 'cigar', 'ref_name', 'ref_start', 'ref_end', 
        'is_reverse', 'mapq', 'as_score', 'strand', 
        'qs', 'qe', 'density', 'is_valid', 'is_backbone'
    ]
                 
    def __init__(self, read: Any = None, from_dict: Optional[Dict] = None):
        self.q_name = ""
        self.flag = 0
        self.cigar = []  
        self.ref_name = ""
        self.ref_start = 0
        self.ref_end = 0
        self.is_reverse = False
        self.mapq = 0
        self.as_score = 0.0
        self.is_backbone = False 
        
        is_original = False 
        
        if read:
            self._init_from_read(read)
            is_original = True
        elif from_dict:
            self._init_from_dict(from_dict)
            if "qs" not in from_dict:
                is_original = True
            
        if not hasattr(self, 'qs') and is_original:
            self.qs, self.qe, self.strand = self._calc_coords_from_cigar()
            
        # Discount score based on mapping quality
        discount = 1.0 if self.mapq >= 20 else 0.8
        self.as_score *= discount
        
        q_span = self.qe - self.qs
        self.density = (self.as_score / q_span) if q_span > 0 else 0.0
        
        self.is_valid = bool(self.cigar and q_span > 0)

        if is_original and self.is_valid:
            self._refine_edges()

    def _get_as_tag(self, read) -> float:
        """Safely retrieve or estimate the alignment score (AS tag)."""
        try: 
            return float(read.get_tag("AS"))
        except KeyError: 
            return float((read.query_alignment_end - read.query_alignment_start) * MATCH_SCORE_DEFAULT / 2)

    def _init_from_read(self, read):
        """Initialize from a pysam AlignedSegment object."""
        self.q_name = read.query_name or ""
        self.flag = read.flag or 0
        self.cigar = read.cigartuples if read.cigartuples else []
        self.ref_name = read.reference_name or ""
        self.ref_start = read.reference_start or 0
        self.ref_end = read.reference_end or 0
        self.is_reverse = read.is_reverse or False
        self.mapq = read.mapping_quality or 0
        self.as_score = self._get_as_tag(read)
        self.strand = '-' if self.is_reverse else '+'

    def _init_from_dict(self, from_dict: Dict):
        """Initialize from dictionary serialization."""
        self.q_name = from_dict.get("q_name", "")
        self.flag = from_dict.get("flag", 0)
        self.cigar = from_dict.get("cigar", [])
        self.ref_name = from_dict.get("ref_name", "")
        self.ref_start = from_dict.get("ref_start", 0)
        self.ref_end = from_dict.get("ref_end", 0)
        self.is_reverse = from_dict.get("is_reverse", False)
        self.mapq = from_dict.get("mapq", 0)
        self.as_score = float(from_dict.get("as_score", 0.0))
        self.strand = from_dict.get("strand", '-' if self.is_reverse else '+')

        if "qs" in from_dict and "qe" in from_dict:
            self.qs = from_dict["qs"]
            self.qe = from_dict["qe"]
       
    def validate(self):
        """Validates critical alignment fields."""
        assert isinstance(self.cigar, list), f"CIGAR must be a list, current: {type(self.cigar)}"
        assert isinstance(self.ref_start, int), f"ref_start must be an int, current: {type(self.ref_start)}"
        assert self.strand in ['+', '-'], f"strand must be +/- , current: {self.strand}"

    def _calc_coords_from_cigar(self) -> Tuple[int, int, str]:
        """Calculates query start/end mapping coordinates considering strand direction."""
        if not self.cigar: 
            return 0, 0, "+"

        q_ops = {0, 1, 4, 5, 7, 8}
        total_len = sum(l for o, l in self.cigar if o in q_ops)

        left_clip = 0
        for op, length in self.cigar:
            if op in {4, 5}: left_clip += length
            elif op in {0, 1, 7, 8}: break

        aligned_q_len = sum(l for op, l in self.cigar if op in {0, 1, 7, 8})
        bam_qs = left_clip
        bam_qe = left_clip + aligned_q_len

        if self.is_reverse:
            return total_len - bam_qe, total_len - bam_qs, "-"
        return bam_qs, bam_qe, "+"
   
    def _refine_edges(self):
        """Identify stable anchors based on score balance, trimming unstable edges."""
        if not self.cigar: return

        trim_q_left, _ = self._find_trim_length_scoring(self.cigar, reverse_scan=False)
        trim_q_right, _ = self._find_trim_length_scoring(self.cigar, reverse_scan=True)
        
        if trim_q_left == 0 and trim_q_right == 0:
            return

        total_q = sum(l for o, l in self.cigar if o in {0, 1, 4, 5, 7, 8})

        if not self.is_reverse:
            target_qs = trim_q_left
            target_qe = total_q - trim_q_right
        else:
            target_qs = trim_q_right
            target_qe = total_q - trim_q_left

        if target_qs >= target_qe:
            self.is_valid = False
            return
        if target_qs == self.qs and target_qe == self.qe:
            return
            
        refined_seg = self.split_at_query_coords(target_qs, target_qe)
        
        if refined_seg:
            self.cigar = refined_seg.cigar
            self.ref_start = refined_seg.ref_start
            self.ref_end = refined_seg.ref_end
            self.qs = refined_seg.qs
            self.qe = refined_seg.qe
            self.as_score = refined_seg.as_score
            self.density = refined_seg.density
        else:
            self.is_valid = False

    def _find_trim_length_scoring(self, cigar: List[Tuple[int, int]], reverse_scan: bool = True) -> Tuple[int, int]:
        HIGH_ANCHOR = 15      
        safe_index = -1
        
        iterator = range(len(cigar)-1, -1, -1) if reverse_scan else range(len(cigar))
        ops = list(iterator) 
        
        for i in range(len(ops)):
            curr_cigar_idx = ops[i]
            op, length = cigar[curr_cigar_idx]
            
            if op in {0, 7} and length >= MIN_ANCHOR_LEN:
                current_score = length * SCORE_MATCH
                is_safe_anchor = True
                
                if length < HIGH_ANCHOR:
                    for k in range(1, ANCHOR_LOOKAHEAD + 1):
                        if i + k >= len(ops): break 
                        
                        nop, nlen = cigar[ops[i + k]]
                        if nop in {0, 7}:      current_score += nlen * SCORE_MATCH
                        elif nop in {1, 2, 8}: current_score += nlen * SCORE_MISMATCH
                        elif nop in {4, 5}:    current_score = -999 
                        
                    if current_score < 0:
                        is_safe_anchor = False
                
                if is_safe_anchor:
                    safe_index = curr_cigar_idx
                    break 
            
        if safe_index == -1:
            return 0, 0 
            
        removed_q, removed_r = 0, 0
        scan_range = range(safe_index + 1, len(cigar)) if reverse_scan else range(0, safe_index)
            
        for k in scan_range:
            op, length = cigar[k]
            if op in {0, 1, 4, 5, 7, 8}: removed_q += length
            if op in {0, 2, 3, 7, 8}: removed_r += length
            
        return removed_q, removed_r

    def split_at_query_coords(self, split_qs: int, split_qe: int) -> Optional['AlignedSegment']:
        """Slices the alignment segment given a target query coordinate interval."""
        q_ops = {0, 1, 4, 5, 7, 8}
        total_len = sum(l for o, l in self.cigar if o in q_ops)
        
        target_qs_fwd = total_len - split_qe if self.is_reverse else split_qs
        target_qe_fwd = total_len - split_qs if self.is_reverse else split_qe

        # Smart anchor snapping
        snap_qs, snap_qe = target_qs_fwd, target_qe_fwd
        tmp_q = 0
        for op, length in self.cigar:
            if op not in q_ops: continue
            op_end = tmp_q + length
            if op in {1, 4, 5}: 
                if tmp_q < snap_qs < op_end: snap_qs = op_end
                if tmp_q < snap_qe < op_end: snap_qe = tmp_q
            tmp_q += length

        target_qs_fwd, target_qe_fwd = snap_qs, snap_qe
        if target_qs_fwd >= target_qe_fwd: return None 

        new_cigar = []
        new_ref_start = self.ref_start
        curr_q = 0  
        curr_r = self.ref_start
        found_start = False
        kept_q_len = 0 

        for op, length in self.cigar:
            consumes_q = op in {0, 1, 4, 5, 7, 8}
            consumes_r = op in {0, 2, 3, 7, 8}
            
            if not consumes_q:
                if target_qs_fwd <= curr_q < target_qe_fwd:
                    if not found_start:
                        new_ref_start = curr_r
                        found_start = True
                    new_cigar.append((op, length))
                if consumes_r: curr_r += length
                continue

            op_q_start, op_q_end = curr_q, curr_q + length
            overlap_start = max(op_q_start, target_qs_fwd)
            overlap_end = min(op_q_end, target_qe_fwd)
            
            if overlap_end > overlap_start:
                kept_len = overlap_end - overlap_start
                if not found_start:
                    new_ref_start = (curr_r + overlap_start - op_q_start) if consumes_r else curr_r
                    found_start = True
                new_cigar.append((op, kept_len))
                if op in {0, 1, 7, 8}: kept_q_len += kept_len
            
            curr_q += length
            if consumes_r: curr_r += length
            if curr_q >= target_qe_fwd: break
                
        # Edge trimming
        for idx_seq in (range(len(new_cigar)), range(len(new_cigar)-1, -1, -1)):
            for i in idx_seq:
                op, length = new_cigar[i]
                if op == 1: new_cigar[i] = (4, length) 
                elif op in {0, 2, 3, 7, 8}: break      

        def _calc_unaligned(cigar_iter) -> int:
            unaligned = 0
            for op, length in cigar_iter:
                if op in {4, 5}: unaligned += length
                elif op in {0, 1, 7, 8}: break
            return unaligned

        aligned_qs_fwd = min(target_qs_fwd + _calc_unaligned(new_cigar), target_qe_fwd)
        aligned_qe_fwd = max(target_qe_fwd - _calc_unaligned(reversed(new_cigar)), target_qs_fwd)

        if self.is_reverse:
            final_qs, final_qe = total_len - aligned_qe_fwd, total_len - aligned_qs_fwd
        else:
            final_qs, final_qe = aligned_qs_fwd, aligned_qe_fwd
            
        ref_consumed = sum(l for o, l in new_cigar if o in {0, 2, 3, 7, 8})
        
        orig_valid_len = self.qe - self.qs
        new_as_score = float(self.as_score * min(1.0, kept_q_len / orig_valid_len)) if orig_valid_len > 0 and kept_q_len > 0 else 0.0
        
        return AlignedSegment(from_dict={
            "q_name": self.q_name, "flag": self.flag, "cigar": new_cigar,
            "ref_name": self.ref_name, "ref_start": new_ref_start,
            "ref_end": new_ref_start + ref_consumed, "is_reverse": self.is_reverse,
            "mapq": self.mapq, "as_score": new_as_score,
            "qs": final_qs, "qe": final_qe, 
            "strand": "-" if self.is_reverse else "+"
        })
    
    def to_block_tuple(self) -> Tuple:
        return (self.qs, self.qe, self.ref_start, self.ref_end, self.ref_name, self.strand)

    def parse_intra_alignment_sv(self, query_seq: Optional[str] = None, min_size: int = 50, merge_gap_tol: int = 50) -> List[Dict]:
        """Parses INDELs found directly inside the CIGAR string."""
        svs = []
        if not self.cigar: return svs
        
        curr_ref, curr_bam_seq_idx = self.ref_start, 0 
        
        # Determine left soft/hard clips for sequence coordinate mapping
        left_clip = 0
        for op, length in self.cigar:
            if op in {4, 5}: left_clip += length
            elif op in {0, 1, 7, 8}: break

        p_ins, ins_gap = None, 0
        p_del, del_gap = None, 0

        def flush_ins():
            nonlocal p_ins
            if p_ins and p_ins['len'] >= min_size:
                svs.append({
                    'type': 'INS', 'chr': self.ref_name, 'pos': p_ins['pos'],
                    'len': p_ins['len'], 'read': self.q_name,
                    'info': {"SEQ": p_ins['seq'], "SOURCE": "CIGAR_MERGED" if "MERGED" in p_ins.get('flag', "") else "CIGAR"}
                })
            p_ins = None

        def flush_del():
            nonlocal p_del
            if p_del and p_del['len'] >= min_size:
                svs.append({
                    'type': 'DEL', 'chr': self.ref_name, 'pos': p_del['pos'],
                    'len': p_del['len'], 'read': self.q_name,
                    'info': {"SOURCE": "CIGAR_MERGED" if "MERGED" in p_del.get('flag', "") else "CIGAR"}
                })
            p_del = None

        for op, length in self.cigar:
            consumes_seq = op in {0, 1, 4, 5, 7, 8}
            consumes_ref = op in {0, 2, 3, 7, 8}

            if op == 1: # INS
                ins_seq = ""
                if query_seq:
                    # Robust absolute coordinate mapping resolving reverse strand misalignments
                    if not self.is_reverse:
                        g_start = self.qs + (curr_bam_seq_idx - left_clip)
                        g_end = g_start + length
                        if 0 <= g_start < g_end <= len(query_seq):
                            ins_seq = query_seq[g_start:g_end]
                    else:
                        g_end = self.qe - (curr_bam_seq_idx - left_clip)
                        g_start = g_end - length
                        if 0 <= g_start < g_end <= len(query_seq):
                            ins_seq = query_seq[g_start:g_end].translate(RC_TRANS)[::-1]

                # VCF standard: POS is the last mapped reference base before the insertion
                sv_pos = curr_ref

                if p_ins:
                    if ins_gap <= merge_gap_tol:
                        p_ins['len'] += length
                        p_ins['seq'] += ins_seq 
                        p_ins['flag'] = "MERGED"
                        ins_gap = 0
                    else:
                        flush_ins()
                        p_ins = {'pos': sv_pos, 'len': length, 'seq': ins_seq}
                        ins_gap = 0
                else:
                    p_ins = {'pos': sv_pos, 'len': length, 'seq': ins_seq}
                    ins_gap = 0
                if p_del: del_gap += length

            elif op == 2: # DEL
                sv_pos = curr_ref  # 0-based start of deleted region
                
                if p_del:
                    if del_gap <= merge_gap_tol:
                        p_del['len'] += length
                        p_del['flag'] = "MERGED"
                        del_gap = 0
                    else:
                        flush_del()
                        p_del = {'pos': sv_pos, 'len': length}
                        del_gap = 0
                else:
                    p_del = {'pos': sv_pos, 'len': length}
                    del_gap = 0
                if p_ins: ins_gap += length
                
            else: 
                if consumes_ref and p_ins: ins_gap += length
                if consumes_seq and consumes_ref and p_del: del_gap += length

            if p_ins and ins_gap > merge_gap_tol: flush_ins()
            if p_del and del_gap > merge_gap_tol: flush_del()

            if consumes_seq: curr_bam_seq_idx += length
            if consumes_ref: curr_ref += length

        flush_ins()
        flush_del()
        return svs

# ==============================================================================
# SV Resolvers
# ==============================================================================

# ── 2-A: OrphanCluster dataclass ──────────────────────────────────────────────
@dataclass
class OrphanCluster:
    """
    condense_orphan_blocks 的输出单元。

    subtype 命名规则（两层，冒号分隔）：
      正向：DUP:TANDEM | DUP:PROXIMAL | DUP:DISP:INTRA | DUP:DISP:INTER
      反向：Inv 前缀，如 InvDUP:TANDEM | InvDUP:DISP:INTRA

    FIX: 新增 true_q_len 和 sub_spans 两个字段。
      - true_q_len：各子块 (qe - qs) 之和，表示孤儿块在 query 中实际消耗的碱基数，
        用于散在重复 INS 和倒位重复 CSV 的 SVLEN，避免用 max(qe)-min(qs) 高估。
      - sub_spans：各子块在参考基因组上的 (ref_start, ref_end) 列表，
        用于在 INFO 字段中精确记录 source 坐标，防止多块聚合后 span 虚大误导用户。
    """
    ref_name  : str
    strand    : str
    ref_start : int
    ref_end   : int
    q_start   : int
    q_end     : int
    raw_type  : str   # "CSV"（同染色体）或 "BND"（跨染色体）
    subtype   : str   # 机制字符串，见上方规则
    blocks    : list = dc_field(default_factory=list, repr=False)
    true_q_len: int  = 0    #
    sub_spans : list = dc_field(default_factory=list, repr=False)  # [(ref_start, ref_end(qs,qe,strand)), ...]
    cluster_spans:Tuple[str,int, int, int, int, str] = () #(chr,rs,re,qs,qe,strand)

    @property
    def ref_len(self) -> int:
        return self.ref_end - self.ref_start

    @property
    def q_len(self) -> int:
        """Span-based query length (max qe - min qs). May include inter-block gaps."""
        return self.q_end - self.q_start

    @property
    def is_inter_chr(self) -> bool:
        return self.raw_type == "BND"

    @property
    def is_inverted(self) -> bool:
        return self.subtype.startswith("Inv")

    @property
    def is_fragmented(self) -> bool:
        """True when sub-blocks have ref gaps > 100 bp between them (non-contiguous source)."""
        if len(self.sub_spans) <= 1:
            return False
        return any(
            self.sub_spans[k + 1][0] - self.sub_spans[k][1] > 100
            for k in range(len(self.sub_spans) - 1)
        )


def _classify_orphan_cluster(
    c_chr           : str,
    c_strand        : str,
    c_r_start       : int,
    c_r_end         : int,
    anchor_chr      : str,
    anchor_strand   : str,
    anchor_bp_left  : int,
    anchor_bp_right : int,
    tandem_max_dist   : int,
    proximal_max_dist : int,
) -> tuple:
    inv_prefix = "Inv" if (c_strand != anchor_strand) else ""

    if c_chr != anchor_chr:
        return ("BND", f"{inv_prefix}DUP:DISP:INTER")

    gap_left  = min(anchor_bp_left, anchor_bp_right)
    gap_right = max(anchor_bp_left, anchor_bp_right)
    dist_to_gap = max(0, gap_left - c_r_end, c_r_start - gap_right)
    orphan_len  = max(1, c_r_end - c_r_start)

    is_strictly_inside = (c_r_start >= gap_left and c_r_end <= gap_right)

    if is_strictly_inside:
        if c_strand != anchor_strand:
            return ("CSV", "InvDUP")
        else:
            return ("CSV", "LOCAL_RETAIN")

    if dist_to_gap <= tandem_max_dist and dist_to_gap <= orphan_len:
        return ("CSV", f"{inv_prefix}DUP:TANDEM")
    elif dist_to_gap <= proximal_max_dist:
        return ("CSV", f"{inv_prefix}DUP:PROXIMAL")
    else:
        return ("CSV", f"{inv_prefix}DUP:DISP:INTRA")

def split_segments(segments: List[AlignedSegment]) -> List[AlignedSegment]:
    """Handles parallel query overlaps by iteratively splitting conflict segments."""
    segments.sort(key=lambda x: (x.qe - x.qs)) 
    queue, processed = segments[:], []

    while queue:
        curr = queue.pop(0)
        if not curr.is_valid: continue
        
        curr_consumed = False
        new_processed = []
        
        for existing in processed:
            if not existing.is_valid:
                new_processed.append(existing)
                continue
            
            if existing.qs <= curr.qs and existing.qe >= curr.qe:
                if curr.density > existing.density * SPLIT_SCORE_RATIO:
                    if curr.qs - existing.qs > HALF_REFINE_TOL:
                        left = existing.split_at_query_coords(existing.qs, curr.qs)
                        if left and left.as_score > 0: new_processed.append(left)
                    if existing.qe - curr.qe > HALF_REFINE_TOL:
                        right = existing.split_at_query_coords(curr.qe, existing.qe)
                        if right and right.as_score > 0: new_processed.append(right)
                    existing.is_valid = False 
                else:
                    new_processed.append(existing)

            elif curr.qs <= existing.qs and curr.qe >= existing.qe:
                if existing.density > curr.density * SPLIT_SCORE_RATIO:
                    if existing.qs - curr.qs > REFINE_TOL:
                        left = curr.split_at_query_coords(curr.qs, existing.qs)
                        if left and left.as_score > 0: queue.append(left)
                    if curr.qe - existing.qe > REFINE_TOL:
                        right = curr.split_at_query_coords(existing.qe, curr.qe)
                        if right and right.as_score > 0: queue.append(right)
                    curr_consumed = True
                    new_processed.append(existing)
                    break 
                else:
                    new_processed.append(existing)
            else:
                new_processed.append(existing)
        
        processed = new_processed
        if not curr_consumed:
            processed.append(curr)
            
    return [s for s in processed if s.is_valid]

def solve_chain_dp(nodes: List[AlignedSegment], query_length: int, K: int = 32) -> List[AlignedSegment]:
    """Dynamic programming algorithm to find the optimal segment chain along the read."""
    if not nodes: return []
    
    nodes = sorted(nodes, key=lambda x: (x.qs, x.qe))
    n = len(nodes)
    dp_score, parent = [0.0] * n, [-1] * n
    
    effective_scores = [max(n.as_score, (n.qe - n.qs) * MIN_AS_PER_BP) for n in nodes]

    for i in range(n):
        node_i, eff_as_i = nodes[i], effective_scores[i]
        dp_score[i], best_prev, max_chain_score = eff_as_i, -1, eff_as_i
        
        for j in range(max(0, i - K), i):
            node_j = nodes[j]
            current_score = dp_score[j] + eff_as_i
            penalty = 0.0
            
            q_gap = node_i.qs - node_j.qe
            
            if q_gap < 0:
                q_overlap = abs(q_gap)
                density = eff_as_i / (node_i.qe - node_i.qs) if (node_i.qe - node_i.qs) > 0 else 1.0
                penalty += q_overlap * density * OVERLAP_PENALTY_RATE
                if q_overlap > MAX_ALLOWED_OVERLAP:
                    penalty += (q_overlap - MAX_ALLOWED_OVERLAP) * 0.5
            elif q_gap > 0:
                penalty += q_gap * UNMAPPED_RATE * 0.5
                
            if node_i.ref_name != node_j.ref_name:
                penalty += 500.0 if (node_i.qe - node_i.qs) >= 1000 else 600.0
            elif node_i.strand != node_j.strand:
                penalty += PENALTY_INV
            else:
                r_gap = (node_i.ref_start - node_j.ref_end) if node_i.strand == '+' else (node_j.ref_start - node_i.ref_end) 
                dist = abs(q_gap - r_gap)
                penalty += MAX_PENALTY_SV if dist > MAX_LINEAR_GAP else dist * LINEAR_GAP_RATE
            
            score_from_j = current_score - penalty
            if score_from_j > max_chain_score:
                max_chain_score, best_prev = score_from_j, j
                
        dp_score[i], parent[i] = max_chain_score, best_prev

    best_global_score, best_end_node = -float('inf'), -1
    
    for i in range(n):
        curr = i
        while parent[curr] != -1: curr = parent[curr]
        
        global_score = dp_score[i] - (nodes[curr].qs + query_length - nodes[i].qe) * UNMAPPED_RATE
        if global_score > best_global_score:
            best_global_score, best_end_node = global_score, i
            
    if best_end_node == -1: return []
        
    chain, curr = [], best_end_node
    while curr != -1:
        chain.append(nodes[curr])
        curr = parent[curr]
        
    return chain[::-1]

def find_collinear_backbone(segs: List[AlignedSegment]):
    """Strict Collinear Backbone: Identifies the dominant, uni-directional mapping chain."""
    n = len(segs)
    if n <= 1:
        if n == 1: segs[0].is_backbone = True
        return set(range(n)), segs[0].strand if n == 1 else '+'
        
    # 按照 Query 坐标排序，保证 DP 顺序的绝对正确
    indexed_segs = sorted(enumerate(segs), key=lambda x: (x[1].qs, x[1].qe))
    
    # 投票决定全局主干方向（只允许一个方向成为主干）
    sum_plus = sum((s.qe - s.qs) for s in segs if s.strand == '+')
    sum_minus = sum((s.qe - s.qs) for s in segs if s.strand == '-')
    global_strand = '+' if sum_plus >= sum_minus else '-'
    
    dp, prev_node = [0.0] * n, [-1] * n
    MAX_TOLERABLE_OVERLAP = 1000 
    
    for i in range(n):
        orig_idx_i, seg_i = indexed_segs[i]
        
        # 【核心约束 1】：方向不一致的片段，直接剥夺主干资格，强行沦为游离块 (Orphan)
        if seg_i.strand != global_strand:
            dp[i] = 0.0
            continue
            
        q_len_i = seg_i.qe - seg_i.qs
        dp[i] = q_len_i
        
        for j in range(i):
            orig_idx_j, seg_j = indexed_segs[j]
            if dp[j] <= 0.0: continue
            
            # 【核心约束 2】：跨染色体的片段绝不能连在同一个主干上
            if seg_i.ref_name != seg_j.ref_name:
                continue
                
            ref_overlap = max(0, min(seg_i.ref_end, seg_j.ref_end) - max(seg_i.ref_start, seg_j.ref_start))
            q_overlap = max(0, seg_j.qe - seg_i.qs)
            
            # 严格几何方向的距离计算（因为已知均为 global_strand）
            r_gap = (seg_i.ref_start - seg_j.ref_end) if global_strand == '+' else (seg_j.ref_start - seg_i.ref_end)
            # =================================================================
            # 【自适应重叠率过滤】：替换原本单一的 MAX_TOLERABLE_OVERLAP = 1000
            # =================================================================
            len_i_q = seg_i.qe - seg_i.qs
            len_j_q = seg_j.qe - seg_j.qs
            min_q_len = min(len_i_q, len_j_q)
            
            len_i_r = seg_i.ref_end - seg_i.ref_start
            len_j_r = seg_j.ref_end - seg_j.ref_start
            min_r_len = min(len_i_r, len_j_r)
            
            q_overlap_ratio = q_overlap / min_q_len if min_q_len > 0 else 0
            ref_overlap_ratio = ref_overlap / min_r_len if min_r_len > 0 else 0
            
            MIN_OVERLAP_THRES   = getattr(SVConfig, "BB_OVERLAP_MIN_BASE", 200)
            MAX_OVERLAP_CEILING = getattr(SVConfig, "BB_OVERLAP_MAX_CEILING", 5000)
            OVERLAP_RATIO_CUT   = getattr(SVConfig, "BB_OVERLAP_RATIO_THRES", 0.5)
            
            # 拦截条件判定
            is_q_invalid = (q_overlap > MAX_OVERLAP_CEILING) or (q_overlap > MIN_OVERLAP_THRES and q_overlap_ratio > OVERLAP_RATIO_CUT)
            is_r_invalid = (ref_overlap > MAX_OVERLAP_CEILING) or (ref_overlap > MIN_OVERLAP_THRES and ref_overlap_ratio > OVERLAP_RATIO_CUT)
            is_gap_invalid = (r_gap < -MAX_OVERLAP_CEILING) or (r_gap < -MIN_OVERLAP_THRES and (abs(r_gap) / min_r_len) > OVERLAP_RATIO_CUT)
            
            if is_q_invalid or is_r_invalid or is_gap_invalid:
                continue 
            
            # 去除冗余重叠长度，防止得分虚高
            score = dp[j] + q_len_i - max(ref_overlap, q_overlap)
            if r_gap > 0: score -= min(1000.0, r_gap * 0.001)

            if score > dp[i]:
                dp[i], prev_node[i] = score, j
                
    max_i = max(range(n), key=lambda x: dp[x])
    backbone_indices = set()
    
    curr = max_i
    while curr != -1:
        orig_idx, seg = indexed_segs[curr]
        seg.is_backbone = True 
        backbone_indices.add(orig_idx)
        curr = prev_node[curr]
        
    return backbone_indices, global_strand

def condense_orphan_blocks(
    orphans         : List[AlignedSegment],
    anchor_chr      : str,
    anchor_strand   : str,
    anchor_bp_left  : int,
    anchor_bp_right : int,
    max_ref_dist    : int = 100_000,
    max_query_dist  : int = 50_000,
) -> List[OrphanCluster]:
    """
    将孤儿比对块聚合为空间连贯的簇，然后对每个簇分类。

    聚合条件（全部满足）：
      - 同一参考染色体
      - 同一链方向（不同链代表不同机制，不能合并）
      - ref 边缘距离 ≤ max_ref_dist
      - query 连续性：query 间隔 ≤ max_query_dist
    """
    if not orphans:
        return []

    TANDEM_MAX   = getattr(SVConfig, "TANDEM_MAX_REF_DIST",   500)
    PROXIMAL_MAX = getattr(SVConfig, "PROXIMAL_MAX_REF_DIST", 10_000)

    def _finalize(c_chr, c_strand, c_r_start, c_r_end, c_q_start, c_q_end, c_blocks):
        raw_type, subtype = _classify_orphan_cluster(
            c_chr, c_strand, c_r_start, c_r_end,
            anchor_chr, anchor_strand,
            anchor_bp_left, anchor_bp_right,
            TANDEM_MAX, PROXIMAL_MAX,
        )
        # cre= max(b.ref_end for b in c_blocks)
        # crs=min(b.ref_start for b in c_blocks)
        # cqe= max(b.qe for b in c_blocks)
        # cqs=min(b.qs for b in c_blocks)
        # true_q_len =cqe-cqs
        true_q_len=c_q_end-c_q_start
        # FIX: 记录各子块的 ref 坐标，防止多块合并后 span 虚大
        sub_spans  = [(b.ref_start, b.ref_end,b.qs,b.qe,b.strand) for b in c_blocks]
        # cluster_spans=(c_chr,crs,cre,cqs,cqe,c_strand)
        cluster_spans=(c_chr,c_r_start,c_r_end,c_q_start,c_q_end,c_strand)
        return OrphanCluster(
            ref_name=c_chr, strand=c_strand,
            ref_start=c_r_start, ref_end=c_r_end,
            q_start=c_q_start, q_end=c_q_end,
            raw_type=raw_type, subtype=subtype,
            blocks=c_blocks,
            true_q_len=true_q_len,
            sub_spans=sub_spans,
            cluster_spans=cluster_spans
        )

    clusters = []
    seg0 = orphans[0]
    curr_chr, curr_strand       = seg0.ref_name, seg0.strand
    curr_r_start, curr_r_end    = seg0.ref_start, seg0.ref_end
    curr_q_start, curr_q_end    = seg0.qs, seg0.qe
    curr_blocks                 = [seg0]

    for seg in orphans[1:]:
        ref_edge_dist = min(
            abs(seg.ref_start - curr_r_end),
            abs(seg.ref_end   - curr_r_start),
        )
        mergeable = (
            seg.ref_name == curr_chr
            and seg.strand == curr_strand
            and ref_edge_dist         <= max_ref_dist
            and (seg.qs - curr_q_end) <= max_query_dist
        )
        if mergeable:
            curr_r_start = min(curr_r_start, seg.ref_start)
            curr_r_end   = max(curr_r_end,   seg.ref_end)
            curr_q_end   = max(curr_q_end,   seg.qe)
            curr_blocks.append(seg)
        else:
            clusters.append(_finalize(curr_chr, curr_strand,
                                      curr_r_start, curr_r_end,
                                      curr_q_start, curr_q_end, curr_blocks))
            curr_chr, curr_strand    = seg.ref_name, seg.strand
            curr_r_start, curr_r_end = seg.ref_start, seg.ref_end
            curr_q_start, curr_q_end = seg.qs, seg.qe
            curr_blocks              = [seg]

    clusters.append(_finalize(curr_chr, curr_strand,
                              curr_r_start, curr_r_end,
                              curr_q_start, curr_q_end, curr_blocks))
    return clusters


def resolve_chain_svs(chains: List['AlignedSegment'], main_strand: str, main_idx_set: Set[int], read_name: str, query_seq: Optional[str]) -> List['SVCandidate']:
    """Master controller that extracts linear/complex SVs from alignment chain transitions (Strict Backbone)."""
    candidates = []
    complex_cluster_idx = 0
    clean_block_id = 0
    
    if not main_idx_set or len(main_idx_set) < 2: 
        return candidates
    
    main_idx = sorted(list(main_idx_set))
    
    def extract_insertion_seq(qs: int, qe: int, strand: str) -> str:
        if not query_seq: return ""
        seq_fwd = query_seq[qs:qe]
        return seq_fwd.translate(RC_TRANS)[::-1] if strand == '-' else seq_fwd

    # =======================================================================
    # Shared Linear SV Resolver
    # =======================================================================
    def check_and_extract_linear_sv(b1: 'AlignedSegment', b2: 'AlignedSegment',
                                    chrom: str, strand: str,
                                    source_label: str, isbackbone=True, extra_info: dict = None) -> Optional['SVCandidate']:
        if b1.strand != b2.strand or b1.ref_name != b2.ref_name:
            return None

        dist_q = b2.qs - b1.qe
        if strand == '+': ref_break_1, ref_break_2 = b1.ref_end, b2.ref_start
        else:             ref_break_1, ref_break_2 = b1.ref_start, b2.ref_end

        dist_r = (ref_break_2 - ref_break_1) if strand == '+' else (ref_break_1 - ref_break_2)
        left_bp, right_bp = min(ref_break_1, ref_break_2), max(ref_break_1, ref_break_2)
        sv_len = abs(dist_q - dist_r)
        
        base_info = {"SOURCE": source_label}
        if extra_info: base_info.update(extra_info)

        if abs(dist_q) < getattr(SVConfig, 'GAP_TOLERANCE', 50) and dist_r < -getattr(SVConfig, 'MIN_SV_SIZE', 30):
            if "COMPLEX" in source_label:

                ref_s1, ref_e1 = min(b1.ref_start, b1.ref_end), max(b1.ref_start, b1.ref_end)
                ref_s2, ref_e2 = min(b2.ref_start, b2.ref_end), max(b2.ref_start, b2.ref_end)
                # 计算两段比对块在 Reference 上的重叠长度
                ref_overlap = max(0, min(ref_e1, ref_e2) - max(ref_s1, ref_s2))

                # 如果两段序列在参考基因组上【完全没有重叠】，且往回跳的距离极大（
                # 物理上不可能是一个串联重复，而是倒位断点两侧的序列错乱
                if ref_overlap == 0 or abs(dist_r) > max(3*(ref_e2-ref_s2),1000000) :
                    return None
            dup_len = abs(dist_r) + dist_q 
            filter_tag = "OVERSIZE" if (getattr(SVConfig, 'MAX_SV_SIZE', -1) != -1 and dup_len >= getattr(SVConfig, 'MAX_SV_SIZE', -1)) else None
            sv_type = "INS" if dup_len < getattr(SVConfig, "TANDEM_DUP_INS_THRESHOLD", 1_000) else "DUP"
            pos = right_bp if sv_type == "INS" else left_bp
            if sv_type == "INS":
                base_info["TANDEM_SOURCE"] = f"{chrom}:{left_bp}-{right_bp}"
                
                # ==========================================================
                # 💡 新增：提取短串联重复的插入序列 SEQ
                # ==========================================================
                # 防止由于超长序列导致内存溢出，加入和标准 INS 一样的长度限制保护
                if query_seq :
                    ins_qs = min(b1.qe, b2.qs)
                    # 确保不越界
                    ins_qe = max(0, min(ins_qs + dup_len, len(query_seq)))
                    if ins_qe > ins_qs:
                        base_info["SEQ"] = extract_insertion_seq(ins_qs, ins_qe, strand)
                # ==========================================================
                
            else:
                base_info["END"] = right_bp

            base_info["SUBTYPE"] = "DUP:TANDEM"
            cand = SVCandidate(sv_type, chrom, pos, dup_len, read_name, info=base_info, filter_str=filter_tag)
            cand.q_start, cand.q_end = min(b1.qe, b2.qs), min(b1.qe, b2.qs) + dup_len
            return cand
            # pos = right_bp if sv_type == "INS" else left_bp
            # if sv_type == "INS":
            #     base_info["TANDEM_SOURCE"] = f"{chrom}:{left_bp}-{right_bp}"
            # else:
            #     base_info["END"] = right_bp

            # base_info["SUBTYPE"] = "DUP:TANDEM"
            # cand = SVCandidate(sv_type, chrom, pos, dup_len, read_name, info=base_info, filter_str=filter_tag)
            # cand.q_start, cand.q_end = min(b1.qe, b2.qs), min(b1.qe, b2.qs) + dup_len
            # return cand

        if sv_len < 20: 
            return None

        if sv_len >= getattr(SVConfig, 'MIN_SV_SIZE', 30) + getattr(SVConfig, 'INSDEL_TOLERANCE', 10):
            if abs(dist_q) < getattr(SVConfig, 'GAP_TOLERANCE', 50) or abs(dist_r) < getattr(SVConfig, 'GAP_TOLERANCE', 50):
                if abs(dist_r) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000): return None
                sv_type = "DEL" if dist_r > dist_q else "INS"
                
                if sv_type == "INS" and (getattr(SVConfig, 'ENABLE_LARGE_SEQUENCE', False) or sv_len < getattr(SVConfig, 'MAX_ALLELE_LEN', 5000)):
                    if query_seq:
                        ins_qs, ins_qe = min(b1.qe, b2.qs), max(0, min(min(b1.qe, b2.qs) + sv_len, len(query_seq)))
                        base_info["SEQ"] = extract_insertion_seq(ins_qs, ins_qe, strand)
                
                cand = SVCandidate(sv_type, chrom, left_bp, sv_len, read_name, info=base_info)
                cand.q_start, cand.q_end = min(b1.qe, b2.qs), max(b1.qe, b2.qs)
                return cand
           
            else:
                if abs(dist_r) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000): return None
                if isbackbone and dist_q > 0 and dist_r > 0:
                    ins_len, del_len = int(dist_q), int(dist_r)
                    if ins_len >= getattr(SVConfig, 'GAP_TOLERANCE', 50) and del_len >= getattr(SVConfig, 'GAP_TOLERANCE', 50):
                        if sv_len >= getattr(SVConfig, 'MIN_CSV_SIZE', getattr(SVConfig, 'MIN_SV_SIZE', 30) * 1.2):
                            base_info.update({"SUBTYPE": "INS_DEL", "INS_LEN": ins_len, "DEL_LEN": del_len, "END": left_bp + del_len})
                            if query_seq:
                                ins_qs, ins_qe = min(b1.qe, b2.qs), max(0, min(min(b1.qe, b2.qs) + ins_len, len(query_seq)))
                                base_info["SEQ"] = extract_insertion_seq(ins_qs, ins_qe, strand)
                            
                            cand = SVCandidate("CSV", chrom, left_bp, sv_len, read_name, info=base_info)
                            cand.q_start, cand.q_end = min(b1.qe, b2.qs), max(b1.qe, b2.qs)
                            return cand
        return None

    # =======================================================================
    # Shared Helper: Extract & Merge Nested Linear SVs
    # =======================================================================
    def _extract_and_merge_nested_nodes(nodes: List['AlignedSegment'], source_label: str, macro_event_id: str, nested_dict: dict):
        """提取内部嵌套微变异，统一绑定宏观 EVENT ID，并动态分配唯一 Candidate ID"""
        for k in range(len(nodes) - 1):
            n1, n2 = nodes[k], nodes[k+1]
            if n1.ref_name != n2.ref_name: 
                continue 
            
            # 只注入全局宏观 EVENT
            internal_cand = check_and_extract_linear_sv(
                b1=n1, b2=n2, chrom=n1.ref_name, strand=n1.strand, 
                source_label=source_label, isbackbone=False, 
                extra_info={"EVENT": macro_event_id}
            )
            if not internal_cand:
                continue

            merge_key = next((_mk for _mk in nested_dict if internal_cand.raw_type == _mk[0] and abs(internal_cand.pos - _mk[1]) <= 10), None)
            
            if merge_key is None:
                nest_idx = len(nested_dict)
                # 为该嵌套变异分配唯一的 VCF ID
                internal_cand.id = f"{macro_event_id}_NEST_{nest_idx}"
                internal_cand.info["MERGED_BLOCKS"] = 1
                nested_dict[(internal_cand.raw_type, internal_cand.pos)] = internal_cand
            else:
                existing = nested_dict[merge_key]
                e_span = existing.q_end - existing.q_start
                n_span = internal_cand.q_end - internal_cand.q_start
                
                is_same_detection = True
                if e_span > 0 and n_span > 0:
                    overlap = max(0, min(existing.q_end, internal_cand.q_end) - max(existing.q_start, internal_cand.q_start))
                    shorter = max(1, min(e_span, n_span))
                    is_same_detection = (overlap / shorter) > 0.5

                if is_same_detection:
                    if internal_cand.sv_len > existing.sv_len:
                        existing.sv_len = internal_cand.sv_len
                        for key in ["INS_LEN", "DEL_LEN", "SEQ", "END"]:
                            if key in internal_cand.info: existing.info[key] = internal_cand.info[key]
                        existing.end = internal_cand.end
                else:
                    existing.sv_len += internal_cand.sv_len
                    existing.info["MERGED_BLOCKS"] = existing.info.get("MERGED_BLOCKS", 1) + 1
                    existing.info["INS_LEN"] = existing.info.get("INS_LEN", 0) + internal_cand.info.get("INS_LEN", 0)
                    existing.info["DEL_LEN"] = existing.info.get("DEL_LEN", 0) + internal_cand.info.get("DEL_LEN", 0)
                    
                    if internal_cand.info.get("SUBTYPE") == "DUP:TANDEM":
                        existing.info["SUBTYPE"] = "MULTI_DUP:TANDEM"
                    if "SEQ" in existing.info and "SEQ" in internal_cand.info:
                        existing.info["SEQ"] += internal_cand.info["SEQ"]
                    if "END" in internal_cand.info and "END" in existing.info:
                        existing.info["END"] = max(existing.info["END"], internal_cand.info["END"])
                    existing.q_end = max(existing.q_end, internal_cand.q_end)
                    
                    if "END" in existing.info:
                        existing.end = int(existing.info["END"])
                    elif existing.type in ("BND", "TRA"):
                        existing.end = existing.pos
                    elif existing.type == "INS":
                        existing.end = existing.pos + int(existing.info.get("REPLACED_REF_LEN", 0))
                    else:
                        existing.end = existing.pos + existing.sv_len

    # =======================================================================
    # Edge Orphan Handling
    # ======================================================================= 
    def process_edge_orphans(bb_node: 'AlignedSegment', orphans: List['AlignedSegment'], is_head: bool):
        if not orphans: return
        
        bb_pos = bb_node.ref_start if bb_node.strand == '+' else bb_node.ref_end
        if not is_head:
            bb_pos = bb_node.ref_end if bb_node.strand == '+' else bb_node.ref_start
            
        anchor_chr, anchor_strand = bb_node.ref_name, bb_node.strand
        global_event_id = f"{read_name}_{'HEAD' if is_head else 'TAIL'}_EDGE"
        
        merged_clusters = condense_orphan_blocks(
            orphans, anchor_chr=anchor_chr, anchor_strand=anchor_strand,
            anchor_bp_left=bb_pos, anchor_bp_right=bb_pos
        )
        
        nested_merged_dict = {}
        _extract_and_merge_nested_nodes(orphans, "COMPLEX_EDGE", global_event_id, nested_merged_dict)
        candidates.extend(nested_merged_dict.values())
        
        # 严格使用最大外包络物理跨度
        total_ins_len = max(o.qe for o in orphans) - min(o.qs for o in orphans)
        
        sub_type_parts = []
        for idx, m in enumerate(merged_clusters):
            if len(merged_clusters) == 1 and m.strand != anchor_strand and m.ref_name == anchor_chr:
                inv_span_len = m.ref_end - m.ref_start
                TRA_GAP_MIN = getattr(SVConfig, 'MIN_TRA_SIZE', 50000)
                TRA_RATIO_THRESHOLD = getattr(SVConfig, 'TRA_RATIO_THRESHOLD', 5.0)
                raw_gap = abs(m.ref_start - bb_pos) if is_head else abs(bb_pos - m.ref_end)
                
                if raw_gap > TRA_GAP_MIN or raw_gap > inv_span_len * TRA_RATIO_THRESHOLD:
                    m.subtype = m.subtype.replace("TANDEM", "DISP:INTRA")
                else:
                    m.subtype = "INV"
            
            sub_type_parts.append(m.subtype)
            span_info = {
                "EVENT": global_event_id, # "CPLX_IDX": idx,
                "SRCREG": f"{m.ref_name}:{m.ref_start}-{m.ref_end}",
                "SOURCE": "COMPLEX_EDGE", "SUBTYPE": m.subtype,
                "DEL_LEN": 0, "INS_LEN": m.true_q_len, "END": bb_pos
            }
            if m.is_fragmented:
                span_info["SUBREGS"] = ",".join(f"{s}-{e}" for s, e,qs,qe,strd in m.sub_spans)
            
            span_cand = SVCandidate("CSV", anchor_chr, bb_pos, m.true_q_len, read_name, info=span_info, filter_str="EDGE_COMPLEX")
            span_cand.id = f"{global_event_id}_SPAN_{idx}"
            candidates.append(span_cand)
            
        eventtype = merged_clusters[0].subtype if len(merged_clusters) == 1 else "COMPLEX_" + "_".join(sub_type_parts)
        
        summary_cand = SVCandidate(
            "CSV", anchor_chr, bb_pos, total_ins_len, read_name, 
            info={"EVENT": global_event_id, "EVENTTYPE": eventtype, "SUBTYPE":"COMPLEX", "DEL_LEN": 0, "INS_LEN": total_ins_len, "ORPHAN_CLUSTERS": len(merged_clusters)}, 
            filter_str="EDGE_COMPLEX"
        )
        summary_cand.id = f"{global_event_id}_SUMMARY"
        candidates.append(summary_cand)

        closest_orphan = orphans[-1] if is_head else orphans[0]
        orphan_pos = closest_orphan.ref_end if closest_orphan.strand == '+' else closest_orphan.ref_start
        if not is_head:
            orphan_pos = closest_orphan.ref_start if closest_orphan.strand == '+' else closest_orphan.ref_end
            
        strands = (closest_orphan.strand + anchor_strand) if is_head else (anchor_strand + closest_orphan.strand)
        tra_type = "INTER_CHR" if anchor_chr != closest_orphan.ref_name else "INTRA_CHR"
        
        c1 = SVCandidate("BND", anchor_chr, bb_pos, 0, read_name, info={
            "EVENT": global_event_id, "CHR2": closest_orphan.ref_name, "POS2": orphan_pos, "TRA_TYPE": tra_type, "STRANDS": strands, "SOURCE": "COMPLEX_EDGE"
        },filter_str="EDGE_COMPLEX")
        c2 = SVCandidate("BND", closest_orphan.ref_name, orphan_pos, 0, read_name, info={
            "EVENT": global_event_id, "CHR2": anchor_chr, "POS2": bb_pos, "TRA_TYPE": tra_type, "STRANDS": strands[::-1], "SOURCE": "COMPLEX_EDGE"
        },filter_str="EDGE_COMPLEX")
        
        # 交叉互锁 MATEID
        c1.id = f"{global_event_id}_BND_EDGE_B"
        c2.id = f"{global_event_id}_BND_EDGE_A"
        c1.mateid = c2.id
        c2.mateid = c1.id
        c1.info["MATEID"] = c1.mateid
        c2.info["MATEID"] = c2.mateid
        candidates.extend([c1, c2])

    if main_idx[0] > 0:
        process_edge_orphans(bb_node=chains[main_idx[0]], orphans=chains[:main_idx[0]], is_head=True)
        
    if main_idx[-1] < len(chains) - 1:
        process_edge_orphans(bb_node=chains[main_idx[-1]], orphans=chains[main_idx[-1] + 1:], is_head=False)

    # =======================================================================
    # Main Chain Iteration
    # =======================================================================
    for i in range(len(main_idx) - 1):
        idx_L, idx_R = main_idx[i], main_idx[i+1]
        bb_left, bb_right = chains[idx_L], chains[idx_R]
        
        if not (bb_left.is_backbone and bb_right.is_backbone): continue 
        
        orphans_in_between = chains[idx_L + 1 : idx_R]
        anchor_chr, anchor_strand = bb_left.ref_name, bb_left.strand
        anchor_bp_left  = bb_left.ref_end   if anchor_strand == '+' else bb_left.ref_start
        anchor_bp_right = bb_right.ref_start if anchor_strand == '+' else bb_right.ref_end
        anchor_pos, anchor_end = anchor_bp_left, anchor_bp_right
        dist_q = bb_right.qs - bb_left.qe
        dist_r = (bb_right.ref_start - bb_left.ref_end) if anchor_strand == '+' else (bb_left.ref_start - bb_right.ref_end)

        if orphans_in_between:
            clean_block_id += 1
            merged_clusters = condense_orphan_blocks(
                orphans_in_between, anchor_chr=anchor_chr, anchor_strand=anchor_strand,
                anchor_bp_left=anchor_bp_left, anchor_bp_right=anchor_bp_right,
            )
            
            complex_cluster_idx += 1 
            global_event_id = f"{read_name}_CPLX_{complex_cluster_idx}"
            if len(merged_clusters) == 1:
                m = merged_clusters[0]
                is_inv = (m.strand != anchor_strand)
                is_same_chr = (m.ref_name == anchor_chr)

                if is_inv and is_same_chr:
                    inv_start, inv_end = m.ref_start, m.ref_end
                    inv_span_len = inv_end - inv_start

                    if anchor_strand == '+':
                        pos_outer_left, pos_outer_right = bb_left.ref_end, bb_right.ref_start
                        q_raw_left, q_raw_right = m.q_start - bb_left.qe, bb_right.qs - m.q_end
                    else:
                        pos_outer_left, pos_outer_right = bb_right.ref_end, bb_left.ref_start
                        q_raw_left, q_raw_right = bb_right.qs - m.q_end, m.q_start - bb_left.qe

                    raw_gap_left, raw_gap_right = inv_start - pos_outer_left, pos_outer_right - inv_end
                    delta_left, delta_right = q_raw_left - raw_gap_left, q_raw_right - raw_gap_right
                    true_ins_left,  true_del_left  = max(0,  delta_left),  max(0, -delta_left)
                    true_ins_right, true_del_right = max(0,  delta_right), max(0, -delta_right)
                    total_span = pos_outer_right - pos_outer_left
                    
                    TRA_GAP_MIN = getattr(SVConfig, 'MIN_TRA_SIZE', 50000)
                    TRA_RATIO_THRESHOLD = getattr(SVConfig, 'TRA_RATIO_THRESHOLD', 5.0)
                    
                    is_remote_inv = (
                        abs(raw_gap_left)  > TRA_GAP_MIN or abs(raw_gap_right) > TRA_GAP_MIN or
                        abs(raw_gap_left)  > inv_span_len * TRA_RATIO_THRESHOLD or
                        abs(raw_gap_right) > inv_span_len * TRA_RATIO_THRESHOLD
                    )

                    if is_remote_inv:
                        m.subtype = m.subtype.replace("TANDEM", "DISP:INTRA")
                    else:
                        MIN_SV = getattr(SVConfig, 'MIN_SV_SIZE', 30)

                        # 1. 计算两侧附加变异的总和
                        total_ins = true_ins_left + true_ins_right
                        total_del = true_del_left + true_del_right

                        # 2. 基于总和进行整体定性
                        is_ins = total_ins > MIN_SV
                        is_del = total_del > MIN_SV

                        if   is_ins and is_del: inv_subtype = "INV-INS-DEL"
                        elif is_ins:            inv_subtype = "INV-INS"
                        elif is_del:            inv_subtype = "INV-DEL"
                        else:                   inv_subtype = "INV"

                        # 3. 构建基础 info 字典 (确保在 inv_subtype 计算之后)
                        info = {
                            "EVENT":       global_event_id,
                            "SUBTYPE":     inv_subtype,
                            "INV_SPAN":    int(inv_span_len),
                            "BP_DETAILS":  f"LC:{int(q_raw_left)},RC:{int(q_raw_right)},LF:{int(raw_gap_left)},RF:{int(raw_gap_right)}"
                        }
                        
                        # 4. 动态追加总和字段 (仅在发生相应变异时输出，避免冗余)
                        if is_ins: info["BP_INS"] = int(delta_left+delta_right)
                        if is_del: info["BP_DEL"] = int(-(delta_right+delta_left))
                        
                        
                        if inv_subtype == "INV":
                            inv_cand = SVCandidate("INV", anchor_chr, inv_start, inv_span_len, read_name, info=info)
                            inv_cand.id = f"{global_event_id}_INV"
                            candidates.append(inv_cand)
                        else:
                            macro_ins_len = int(q_raw_left + inv_span_len + q_raw_right)
                            info.update({"DEL_LEN": int(total_span), "INS_LEN": macro_ins_len, "END": int(pos_outer_right)})
                            csv_cand = SVCandidate("CSV", anchor_chr, pos_outer_left, macro_ins_len, read_name, info=info)
                            csv_cand.id = f"{global_event_id}_CSV"
                            candidates.append(csv_cand)
                        
                        if len(m.blocks) > 1:
                            nested_merged_dict = {}
                            _extract_and_merge_nested_nodes(m.blocks, "COMPLEX", global_event_id, nested_merged_dict)
                            candidates.extend(nested_merged_dict.values())
                        continue
                        
            nested_merged_dict = {}
            all_nodes = [bb_left] + orphans_in_between + [bb_right]
            _extract_and_merge_nested_nodes(all_nodes, "COMPLEX", global_event_id, nested_merged_dict)
            candidates.extend(nested_merged_dict.values())
            sub_type_info = []
            
            # 1. 记录 Backbone 的跨度（修正 NameError：直接使用 bb_left.qe 和 bb_right.qs）
            sub_type_info.append(
                f"backbone,{bb_left.ref_name},{anchor_bp_left},{anchor_bp_right},{bb_left.qe},{bb_right.qs},{main_strand}"
            )
            sub_type_parts = []
            prev_chr, prev_strand, prev_pos = anchor_chr, anchor_strand, anchor_pos

            for idx, m in enumerate(merged_clusters):
                sub_type_parts.append(m.subtype)
                sub_type_info.append(
                    f"{m.subtype},{m.cluster_spans[0]},{m.cluster_spans[1]},{m.cluster_spans[2]},{m.cluster_spans[3]},{m.cluster_spans[4]},{m.cluster_spans[5]}"
                )
                m_enter_pos = m.ref_start if m.strand == '+' else m.ref_end
                m_leave_pos = m.ref_end   if m.strand == '+' else m.ref_start

                is_inter_chr     = (prev_chr != m.ref_name)
                is_distant_intra = (not is_inter_chr and abs(prev_pos - m_enter_pos) > getattr(SVConfig, 'MIN_TRA_SIZE', 50_000))
                is_distant_jump = is_inter_chr or is_distant_intra

                bnd_entry = SVCandidate("BND", prev_chr, prev_pos, 0, read_name, info={
                    "EVENT": global_event_id, "CHR2": m.ref_name, "POS2": m_enter_pos,
                    "TRA_TYPE": "INTER_CHR" if is_inter_chr else "INTRA_CHR", "STRANDS": prev_strand + m.strand, "SOURCE": "COMPLEX",
                }, filter_str=None if is_distant_jump else "LOCAL_BND")
                bnd_entry.id = f"{global_event_id}_BND_{idx}_A"
                bnd_entry.mateid = f"{global_event_id}_BND_{idx}_B"

                bnd_exit = SVCandidate("BND", m.ref_name, m_enter_pos, 0, read_name, info={
                    "EVENT": global_event_id, "CHR2": prev_chr, "POS2": prev_pos,
                    "TRA_TYPE": "INTER_CHR" if is_inter_chr else "INTRA_CHR", "STRANDS": m.strand + prev_strand, "SOURCE": "COMPLEX",
                }, filter_str=None if is_distant_jump else "LOCAL_BND")
                bnd_exit.id = bnd_entry.mateid
                bnd_exit.mateid = bnd_entry.id

                bnd_entry.info["MATEID"] = bnd_entry.mateid
                bnd_exit.info["MATEID"] = bnd_exit.mateid

                candidates.extend([bnd_entry, bnd_exit])
                prev_chr, prev_strand, prev_pos = m.ref_name, m.strand, m_leave_pos
                
                # ── 基础 Span 记录构建（保留全部溯源标签） ──
                span_info = {
                    "EVENT": global_event_id, 
                   # "CPLX_IDX": idx,
                    "SRCREG": f"{m.ref_name}:{m.ref_start}-{m.ref_end}",
                    "SUBTYPE": m.subtype, 
                    "DEL_LEN": 0, 
                    "INS_LEN": m.true_q_len
                }
                if m.is_fragmented:
                    span_info["SUBREGS"] = ",".join(f"{s}-{e}" for s, e,qs,qe,strd in m.sub_spans)
                    
                if m.subtype != "DUP:TANDEM":
                    # 1. 局部保留 (LOCAL_RETAIN) 局部混乱
                    if m.subtype == "LOCAL_RETAIN":
                        # m.subtype="DUP"
                        # span_info["END"] = m.ref_end
                        # dup_len = m.ref_end - m.ref_start
                        # span_cand = SVCandidate("CSV", anchor_chr, m.ref_start, dup_len , read_name, info=span_info, filter_str="LOCAL_RETAIN")
                        # span_cand.id = f"{global_event_id}_SPAN_{idx}"
                        # candidates.append(span_cand)
                        pass
                        
                    # 2. 嵌套倒位 (INV) 
                    # 必须使用其真实的参考基因组坐标 (m.ref_start) 和真实倒位跨度
                    elif m.subtype == "InvDUP":
                        span_info["END"] = m.ref_end  # 倒位必须有明确的 END = POS + SVLEN
                        inv_len = m.ref_end - m.ref_start
                        # 坐标使用 m.ref_start，SVLEN 使用 inv_len
                        span_cand = SVCandidate("CSV", m.ref_name, m.ref_start, inv_len, read_name, info=span_info)
                        span_cand.id = f"{global_event_id}_SPAN_{idx}"
                        candidates.append(span_cand)
                        
                    # 3. 散在重复与其他复杂嵌套 (DUP:DISP 等)
                    # 挂载在主干间隙边缘 (anchor_bp_left)
                    else:
                        span_info["END"] = anchor_bp_left
                        span_cand = SVCandidate("CSV", anchor_chr, anchor_bp_left, m.true_q_len, read_name, info=span_info)
                        span_cand.id = f"{global_event_id}_SPAN_{idx}"
                        candidates.append(span_cand)
    
            bb_right_enter = bb_right.ref_start if bb_right.strand == '+' else bb_right.ref_end
            tail_inter     = (prev_chr != bb_right.ref_name)
            tail_distant   = (not tail_inter and abs(prev_pos - bb_right_enter) > getattr(SVConfig, 'MIN_TRA_SIZE', 50_000))
            
            if tail_inter or tail_distant:
                t1 = SVCandidate("BND", prev_chr, prev_pos, 0, read_name, info={
                    "EVENT": global_event_id, "CHR2": bb_right.ref_name, "POS2": bb_right_enter,
                    "TRA_TYPE": "INTER_CHR" if tail_inter else "INTRA_CHR", "STRANDS": prev_strand + bb_right.strand,
                })
                t1.id = f"{global_event_id}_BND_TAIL_A"
                t1.mateid = f"{global_event_id}_BND_TAIL_B"
                
                t2 = SVCandidate("BND", bb_right.ref_name, bb_right_enter, 0, read_name, info={
                    "EVENT": global_event_id, "CHR2": prev_chr, "POS2": prev_pos,
                    "TRA_TYPE": "INTER_CHR" if tail_inter else "INTRA_CHR", "STRANDS": bb_right.strand + prev_strand,
                })
                t2.id = t1.mateid
                t2.mateid = t1.id
                
                t1.info["MATEID"] = t1.mateid
                t2.info["MATEID"] = t2.mateid
                candidates.extend([t1, t2])
            
            if len(merged_clusters) >= 2:
                up_down      = bool(anchor_pos > anchor_end)
                affected_ref = abs(dist_r)
                #strands_str  = (bb_right.strand + anchor_strand) if up_down else (anchor_strand + bb_right.strand)
                eventtype    = (merged_clusters[0].subtype if len(merged_clusters) == 1 else "COMPLEX_" + "_".join(sub_type_parts))
                summary_pos  = anchor_end if up_down else anchor_pos
                summary_end  = summary_pos + affected_ref
                # 3. 将数组合并为 VCF 兼容的字符串格式 (以 | 分隔)
                serialized_sub_svs = "|".join(sub_type_info)
                
                summary_info = {
                    "EVENT": global_event_id, "EVENTTYPE": eventtype, "SUBTYPE": "COMPLEX",
                    "DEL_LEN": int(affected_ref), "INS_LEN": int(max(0, dist_q)), "END": int(summary_end),
                    "ORPHAN_CLUSTERS": len(merged_clusters), #"BACKBONE_STRAND": strands_str, "BACKBONE_FLIP": up_down,
                    "SUB_SVS": serialized_sub_svs # <--- 修正命名并写入字符串
                }

                summary_cand = SVCandidate("CSV", anchor_chr, summary_pos, max(0, dist_q), read_name, info=summary_info, filter_str="COMPLEX")
                summary_cand.id = f"{global_event_id}_SUMMARY"
                candidates.append(summary_cand)
               
        else:
            backbone_cand = check_and_extract_linear_sv(
                b1=bb_left, b2=bb_right, chrom=anchor_chr, strand=anchor_strand, source_label="backbone"
            )
            if backbone_cand: candidates.append(backbone_cand)

    return candidates