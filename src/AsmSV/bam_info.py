# ==============================================================================
# Configuration & Static Constants
# ==============================================================================
from typing import List, Dict, Tuple, Set, Optional, Any
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

                # VCF standard: Position must be the reference base strictly preceding the variant
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
                sv_pos = curr_ref # Base preceding deletion
                
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
    """Identifies the dominant collinear mapping sequence to act as the primary backbone."""
    n = len(segs)
    if n <= 1:
        if n == 1: segs[0].is_backbone = True
        return set(range(n)), segs[0].strand if n == 1 else '+'
        
    sum_plus = sum((s.qe - s.qs) for s in segs if s.strand == '+')
    sum_minus = sum((s.qe - s.qs) for s in segs if s.strand == '-')
    global_strand = '+' if sum_plus >= sum_minus else '-'
    
    dp, prev_node = [0.0] * n, [-1] * n
    MAX_TOLERABLE_OVERLAP = 1000 
    
    for i, seg_i in enumerate(segs):
        q_len_i = seg_i.qe - seg_i.qs
        dp[i] = q_len_i
        
        for j in range(i):
            seg_j = segs[j]
            score = dp[j] + q_len_i
            
            if seg_i.ref_name != seg_j.ref_name:
                score -= (q_len_i + (seg_j.qe - seg_j.qs) + 999999) 
            else:
                ref_overlap = max(0, min(seg_i.ref_end, seg_j.ref_end) - max(seg_i.ref_start, seg_j.ref_start))
                r_gap = (seg_i.ref_start - seg_j.ref_end) if global_strand == '+' else (seg_j.ref_start - seg_i.ref_end)
                
                if ref_overlap > MAX_TOLERABLE_OVERLAP or r_gap < -MAX_TOLERABLE_OVERLAP:
                    score -= (q_len_i + (seg_j.qe - seg_j.qs) + 999999)
                else:
                    if seg_i.strand != seg_j.strand: score -= 500
                    if r_gap > 0: score -= min(1000.0, r_gap * 0.001)

            if score > dp[i]:
                dp[i], prev_node[i] = score, j
                
    max_idx = max(range(n), key=lambda x: dp[x])
    backbone_indices = set()
    
    curr = max_idx
    while curr != -1:
        segs[curr].is_backbone = True 
        backbone_indices.add(curr)
        curr = prev_node[curr]
        
    return backbone_indices, global_strand

def condense_orphan_blocks(orphans: List[AlignedSegment], anchor_chr: str, anchor_pos: int, anchor_strand: str, max_ref_dist: int = 100000, max_query_dist: int = 50000) -> List[Dict]:
    """Condenses complex floating orphan fragments into cluster events."""
    if not orphans: return []
    condensed_clusters = []
    
    curr_chr, curr_strand = orphans[0].ref_name, orphans[0].strand
    curr_r_start, curr_r_end = orphans[0].ref_start, orphans[0].ref_end
    curr_q_start, curr_q_end = orphans[0].qs, orphans[0].qe
    curr_blocks = [orphans[0]] 
    def finalize_cluster(c_chr, c_strand, c_r_start, c_r_end, c_q_start, c_q_end, c_blocks):
        is_inv = (c_strand != anchor_strand)
        prefix = "Inv" if is_inv else ""
        
        if c_chr != anchor_chr: 

            sv_type = "TRA" 
        else:
            dist_to_anchor = min(abs(c_r_start - anchor_pos), abs(c_r_end - anchor_pos))

            if dist_to_anchor > 1000000: 
                sv_type = "DUP:ISP" 
            else:
                is_upstream = (anchor_strand == '+' and c_r_start <= anchor_pos) or \
                              (anchor_strand == '-' and c_r_start >= anchor_pos)
                
                if is_upstream:
                    sv_type = "DUP:TANDEM"
                else:
                    sv_type = "DUP:ISP"
                if is_inv:
                    prefix, sv_type = "", "InvDUP" 

        return {
            "sv_type": f"{prefix}{sv_type}",       
            "ref_name": c_chr, "strand": c_strand,
            "ref_start": c_r_start, "ref_end": c_r_end, "ref_len": c_r_end - c_r_start,
            "q_start": c_q_start, "q_end": c_q_end, "q_len": c_q_end - c_q_start,
            "original_blocks": c_blocks  
        }
    
    for seg in orphans[1:]:
        can_merge = False
        if seg.ref_name == curr_chr and seg.strand == curr_strand:
            ref_dist = min(abs(seg.ref_start - curr_r_end), abs(seg.ref_end - curr_r_start))
            if ref_dist <= max_ref_dist and (seg.qs - curr_q_end) <= max_query_dist:
                can_merge = True
                
        if can_merge:
            curr_r_start, curr_r_end = min(curr_r_start, seg.ref_start), max(curr_r_end, seg.ref_end)
            curr_q_end = max(curr_q_end, seg.qe)
            curr_blocks.append(seg)
        else:
            condensed_clusters.append(finalize_cluster(curr_chr, curr_strand, curr_r_start, curr_r_end, curr_q_start, curr_q_end, curr_blocks))
            curr_chr, curr_strand = seg.ref_name, seg.strand
            curr_r_start, curr_r_end = seg.ref_start, seg.ref_end
            curr_q_start, curr_q_end = seg.qs, seg.qe
            curr_blocks = [seg]
            
    condensed_clusters.append(finalize_cluster(curr_chr, curr_strand, curr_r_start, curr_r_end, curr_q_start, curr_q_end, curr_blocks))
    return condensed_clusters

def pair_inv_signals(signals: List[Dict], read_name: str, main_strand: str) -> List[SVCandidate]:
    """Matches standalone inversion breakpoints into consolidated inversion/translocation events."""
    candidates = []
    if not signals:
        return candidates
        
    paired_indices = set()
    
    TRA_GAP_MIN = SVConfig.MIN_TRA_SIZE
    TRA_RATIO_THRESHOLD = SVConfig.TRA_RATIO_THRESHOLD
    MIN_SV_SIZE = SVConfig.MIN_SV_SIZE

    for i in range(len(signals)):
        if i in paired_indices: continue
        s1 = signals[i]
        
        for j in range(i + 1, len(signals)):
            if j in paired_indices: continue
            s2 = signals[j]
            
            # Topological closure requirement: must be contiguous on the same chromosome
            if s1["chrom"] != s2["chrom"] or s1.get("clean_id", 0) != s2.get("clean_id", -1):
                continue
                
            is_fwd = (s1['type'] == '3to3' and s2['type'] == '5to5')
            is_rev = (s1['type'] == '5to5' and s2['type'] == '3to3')
            if not (is_fwd or is_rev):
                continue
            inv_strand = '-' if is_fwd else '+'
            
            # Ignore false-positive inversions caused by noisy mappings
            if inv_strand == main_strand:
                break
                
            if is_fwd:
                strands = "+-"
                pos_outer_left = s1['pos1']
                inv_end, inv_start = s1['pos2'], s2['pos1']
                pos_outer_right = s2['pos2']
                q_raw_left, q_raw_right = s1['dist_q'], s2['dist_q']
            else:
                strands = "-+"
                pos_outer_left = s2['pos2']
                inv_start, inv_end = s1['pos2'], s2['pos1']
                pos_outer_right = s1['pos1']
                q_raw_left, q_raw_right = s2['dist_q'], s1['dist_q']
                
            start, end = min(inv_start, inv_end), max(inv_start, inv_end)
            inv_span_len = end - start
            total_span = pos_outer_right - pos_outer_left
            
            raw_gap_left = start - pos_outer_left
            raw_gap_right = pos_outer_right - end

            # Delta check: differentiate aligner truncation from real INDEL presence
            delta_left = q_raw_left - raw_gap_left
            delta_right = q_raw_right - raw_gap_right

            true_ins_left = max(0, delta_left)
            true_del_left = max(0, -delta_left)
            true_ins_right = max(0, delta_right)
            true_del_right = max(0, -delta_right)

            sv_type = "INV"
            info = {
                "INV_SPAN": int(inv_span_len),
                "LEFT_INS": int(true_ins_left),
                "RIGHT_INS": int(true_ins_right)
            }
            
            # Scenario 1: Source of INV comes from remote location (Translocation context)
            is_remote_inv = (
                raw_gap_left > TRA_GAP_MIN or 
                raw_gap_right > TRA_GAP_MIN or 
                raw_gap_left > inv_span_len * TRA_RATIO_THRESHOLD or 
                raw_gap_right > inv_span_len * TRA_RATIO_THRESHOLD
            )

            if is_remote_inv:
                sv_type = "BND"
                info.update({
                    "SUBTYPE": "InvTRA",
                    "STRANDS": strands,
                    "FLANK_SPAN": int(total_span)
                })
                
            # Scenario 2: Local DUP-INV (Flanks exhibit extensive regression)
            elif total_span < -MIN_SV_SIZE:
                sv_type = "DUP"
                info.update({"SUBTYPE": "InvDUP"})

            # Scenario 3: Native inversion, checking for accompanying INDELs
            else:
                info.update({
                    "LEFT_DEL": int(true_del_left),
                    "RIGHT_DEL": int(true_del_right)
                })
                is_del = (true_del_left + true_del_right) > MIN_SV_SIZE
                is_ins = (true_ins_left + true_ins_right) > MIN_SV_SIZE

                if is_ins and is_del: subtype = "INV-INS-DEL"
                elif is_ins: subtype = "INV-INS"
                elif is_del: subtype = "INV-DEL"
                else: subtype = "INV"
                info["SUBTYPE"] = subtype

            # Commit the paired variations
            if sv_type == "BND":
                cand = SVCandidate("BND", s1['chrom'], pos_outer_left, 0, read_name, info)
                cand.info["CHR2"] = s1['chrom']
                cand.info["POS2"] = pos_outer_right
                candidates.append(SVCandidate("INV", s1['chrom'], start, inv_span_len, read_name, info))
                candidates.append(cand)
            else:
                candidates.append(SVCandidate(sv_type, s1['chrom'], start, inv_span_len, read_name, info))
            
            paired_indices.update([i, j])
            break  

    # Process unpaired inversion breakends
    for i, s in enumerate(signals):
        if i not in paired_indices:
            strands = "++" if s["type"] == "5to5" else "--"
            candidates.append(SVCandidate("BND", s["chrom"], s["pos1"], 0, read_name, info={
                "SUBTYPE": "UNPAIRED_INV_BREAKEND", 
                "STRANDS": strands,
                "CHR2": s["chrom"], 
                "POS2": s["pos2"],
                "LEFT_INS": int(max(0, s["dist_q"])),
                "RIGHT_INS": int(max(0, s["dist_q"]))
            }))
            
    return candidates

def resolve_chain_svs(chains: List[AlignedSegment], main_strand: str, main_idx_set: Set[int], read_name: str, query_seq: Optional[str]) -> List[SVCandidate]:
    """Master controller that extracts linear/complex SVs from alignment chain transitions."""
    candidates, bp_signals = [], []
    complex_cluster_idx = 0
    clean_block_id = 0
    if not main_idx_set or len(main_idx_set) < 2: return candidates
    main_idx = sorted(list(main_idx_set))

    def extract_insertion_seq(qs: int, qe: int, strand: str) -> str:
        """Helper to obtain properly-oriented insertion sequence from the query string."""
        if not query_seq: return ""
        seq_fwd = query_seq[qs:qe]
        return seq_fwd.translate(RC_TRANS)[::-1] if strand == '-' else seq_fwd

    # =======================================================================
    # Shared Linear SV Resolver (Adhering to DRY principle)
    # Applied to both Backbone transitions and internal Orphan structures.
    # =======================================================================
    def check_and_extract_linear_sv(b1: AlignedSegment, b2: AlignedSegment, 
                                    chrom: str, strand: str, 
                                    source_label: str, extra_info: dict = None) -> Optional[SVCandidate]:
        if b1.strand != b2.strand:
            return None 
            
        dist_q = b2.qs - b1.qe
        
        # Enforce left-aligned VCF coordinate specification
        if strand == '+':
            ref_break_1, ref_break_2 = b1.ref_end, b2.ref_start
        else:
            ref_break_1, ref_break_2 = b1.ref_start, b2.ref_end
            
        dist_r = (ref_break_2 - ref_break_1) if strand == '+' else (ref_break_1 - ref_break_2)
        anchor_pos = min(ref_break_1, ref_break_2)
                
        sv_len = abs(dist_q - dist_r)
        base_info = {"SOURCE": source_label}
        if extra_info: base_info.update(extra_info)
        
        # Check for Tandem Duplications
        if abs(dist_q) < SVConfig.GAP_TOLERANCE and dist_r < -SVConfig.MIN_SV_SIZE:
            dup_len = abs(dist_r) + dist_q
            if SVConfig.MAX_SV_SIZE == -1 or (SVConfig.MAX_SV_SIZE > 0 and dup_len < SVConfig.MAX_SV_SIZE):
               base_info["SUBTYPE"] = "DUP:TANDEM"
            return SVCandidate("INS", chrom, anchor_pos, dup_len, read_name, info=base_info)
       
        if sv_len < 20:
            return None
            
        # Standard INS/DEL processing
        if sv_len >= SVConfig.MIN_SV_SIZE + SVConfig.INSDEL_TOLERANCE:
            if abs(dist_q) < SVConfig.GAP_TOLERANCE or abs(dist_r) < SVConfig.GAP_TOLERANCE:
                # Delegate massive gaps to BND handling
                if abs(dist_r) > SVConfig.MIN_TRA_SIZE:
                    return None
                    
                sv_type = "DEL" if dist_r > dist_q else "INS"
                
                # Retrieve structural insertion sequence logic
                if sv_type == "INS":
                    if sv_len < SVConfig.MAX_ALLELE_LEN or SVConfig.MAX_SV_SIZE == -1:
                        if query_seq:
                            ins_qs = min(b1.qe, b2.qs)
                            ins_qe = max(0, min(ins_qs + sv_len, len(query_seq)))
                            base_info["SEQ"] = extract_insertion_seq(ins_qs, ins_qe, strand)
                                
                return SVCandidate(sv_type, chrom, anchor_pos, sv_len, read_name, info=base_info)
                
        return None

    # ================= Main Chain Iteration =================
    for i in range(len(main_idx) - 1):
        idx_L, idx_R = main_idx[i], main_idx[i+1]
        bb_left, bb_right = chains[idx_L], chains[idx_R]
        
        if not (bb_left.is_backbone and bb_right.is_backbone): continue 
        
        orphans_in_between = chains[idx_L + 1 : idx_R]
        anchor_chr, anchor_strand = bb_left.ref_name, bb_left.strand
        anchor_pos = bb_left.ref_end if anchor_strand == '+' else bb_left.ref_start
        dist_q = bb_right.qs - bb_left.qe
        dist_r = (bb_right.ref_start - bb_left.ref_end) if anchor_strand == '+' else (bb_left.ref_start - bb_right.ref_end)
            
        # =======================================================
        # Scenario A: Presence of orphans (Complex/Jumping Variations)
        # =======================================================
        if orphans_in_between:
            clean_block_id += 1
            merged_clusters = condense_orphan_blocks(orphans_in_between, anchor_chr, anchor_pos, anchor_strand)
            
            complex_cluster_idx += 1 
            event_id = f"{read_name}_CPLX_{complex_cluster_idx}"
            sub_type_str = "COMPLEX"
            
            prev_chr, prev_strand, prev_pos = anchor_chr, anchor_strand, anchor_pos
            
            for idx, m in enumerate(merged_clusters):
                sub_type_str += f"-{m['sv_type']}"
                sub_event_id = f"{event_id}_sub{idx}"
                blocks = m["original_blocks"]

                # Resolve nested events inside the complex cluster
                nested_merged_dict = {}
                for k in range(len(blocks) - 1):
                    internal_cand = check_and_extract_linear_sv(
                        b1=blocks[k], b2=blocks[k+1], 
                        chrom=m["ref_name"], strand=m["strand"], 
                        source_label="COMPLEX_NESTED", 
                        extra_info={"PEVENT": sub_event_id}
                    )
                    
                    if internal_cand:
                        merge_key = None
                        # Tolerance for minor breakpoint wobbling
                        for existing_key in nested_merged_dict.keys():
                            if internal_cand.raw_type == existing_key[0] and abs(internal_cand.pos - existing_key[1]) <= 10:
                                merge_key = existing_key
                                break
                        
                        if merge_key is None:
                            merge_key = (internal_cand.raw_type, internal_cand.pos)
                            internal_cand.info["MERGED_BLOCKS"] = 1
                            nested_merged_dict[merge_key] = internal_cand
                        else:
                            # Concatenate sequence structures of continuous jumps
                            existing_cand = nested_merged_dict[merge_key]
                            existing_cand.sv_len += internal_cand.sv_len
                            
                            if "SEQ" in existing_cand.info and "SEQ" in internal_cand.info:
                                existing_cand.info["SEQ"] += internal_cand.info["SEQ"]
                                
                            if "REPLACED_REF_LEN" in internal_cand.info:
                                existing_cand.info["REPLACED_REF_LEN"] = existing_cand.info.get("REPLACED_REF_LEN", 0) + internal_cand.info["REPLACED_REF_LEN"]
                            
                            existing_cand.info["MERGED_BLOCKS"] = existing_cand.info.get("MERGED_BLOCKS", 1) + 1
                            existing_cand.info["SUBTYPE"] = "MULTIPLE_DUP:TANDEM" if internal_cand.info.get("SUBTYPE") == "DUP:TANDEM" else "MERGED_INS"
                                
                candidates.extend(nested_merged_dict.values())

                m_enter_pos = m["ref_start"] if m["strand"] == '+' else m["ref_end"]
                m_leave_pos = m["ref_end"] if m["strand"] == '+' else m["ref_start"]
                
                is_translocation = (prev_chr != m["ref_name"])
                is_distant_jump = (prev_chr == m["ref_name"] and abs(prev_pos - m_enter_pos) > SVConfig.MIN_TRA_SIZE)
                
                if is_translocation or is_distant_jump:
                    candidates.append(SVCandidate("BND", prev_chr, prev_pos, 0, read_name, info={
                        "PEVENT": event_id, "EVENT": f"{sub_event_id}_BND",
                        "CHR2": m["ref_name"], "POS2": m_enter_pos, 
                        "TRA_TYPE": "INTER_CHR" if is_translocation else "INTRA_CHR",
                        "STRANDS": prev_strand + m["strand"]
                    }))
                
                span_type, sub_info = m["sv_type"], {}
                if span_type == "TRA":
                    span_type, sub_info["SUBTYPE"] = "DUP", "TRA"  
                elif span_type == "InvTRA":
                    span_type, sub_info["SUBTYPE"] = "INV", "InvTRA" 
                sub_info.update({
                    "PEVENT": event_id, "EVENT": sub_event_id,
                    "SOURCE_CHR": m["ref_name"], "SOURCE_POS": m["ref_start"], "SOURCE_END": m["ref_end"],
                    "DEST_CHR": prev_chr, "DEST_POS": prev_pos
                })
                
                candidates.append(SVCandidate(
                    span_type, m["ref_name"], m["ref_start"], m["ref_len"], read_name, info=sub_info
                ))
                
                prev_chr, prev_strand, prev_pos = m["ref_name"], m["strand"], m_leave_pos
                
            bb_right_enter = bb_right.ref_start if bb_right.strand == '+' else bb_right.ref_end
            if prev_chr != bb_right.ref_name or (prev_chr == bb_right.ref_name and abs(prev_pos - bb_right_enter) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000)):
                candidates.append(SVCandidate("BND", prev_chr, prev_pos, 0, read_name, info={
                    "PEVENT": event_id, "EVENT": f"{event_id}_tail_BND",
                    "CHR2": bb_right.ref_name, "POS2": bb_right_enter,
                    "TRA_TYPE": "INTER_CHR" if prev_chr != bb_right.ref_name else "INTRA_CHR",
                    "STRANDS": prev_strand + bb_right.strand
                }))
                
            if len(merged_clusters) >= 2:
                candidates.append(SVCandidate("COMPLEX", anchor_chr, anchor_pos, max(0, dist_q), read_name, info={
                    "EVENT": event_id, "SUBTYPE": sub_type_str,
                    "REPLACED_REF_LEN": int(max(0, dist_r)), "ORPHAN_CLUSTERS": len(merged_clusters),
                    "BACKBONE_STRAND": anchor_strand
                }))

        # =======================================================
        # Scenario B: Clean Backbone Transition (Pure INS/DEL/INV)
        # =======================================================
        else:
            if bb_right.strand != bb_left.strand:
                orient = "3to3" if bb_left.strand == '+' else "5to5"
                exit_pos = bb_left.ref_end if bb_left.strand == '+' else bb_left.ref_start
                enter_pos = bb_right.ref_start if bb_right.strand == '+' else bb_right.ref_end
                
                bp_signals.append({
                    "type": orient, 
                    "chrom": anchor_chr, 
                    "pos1": exit_pos,    
                    "pos2": enter_pos,   
                    "strand_L": bb_left.strand, 
                    "strand_R": bb_right.strand, 
                    "dist_q": dist_q,
                    "clean_id": clean_block_id,
                    "read_name": read_name
                })
                continue
            
            # Extract standard variants utilizing the shared resolving logic
            backbone_cand = check_and_extract_linear_sv(
                b1=bb_left, b2=bb_right, 
                chrom=anchor_chr, strand=anchor_strand, 
                source_label="backbone"
            )
            if backbone_cand:
                candidates.append(backbone_cand)
        
    if bp_signals:
        candidates.extend(pair_inv_signals(bp_signals, read_name, main_strand))
        
    return candidates


