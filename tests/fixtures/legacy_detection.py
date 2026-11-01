"""Conservative legacy split/CSV interpretation; no replacement candidate expansion."""
from typing import List, Dict, Tuple, Set, Optional, Any
from dataclasses import dataclass, field as dc_field
from asmsv_next.config import SVConfig
from asmsv_next.structures import SVCandidate
from asmsv_next.alignment import AlignedSegment
RC_TRANS = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')

@dataclass
class OrphanCluster:
    """Legacy branch summary; true_q_len retains the historical outer-span heuristic."""
    ref_name: str
    strand: str
    ref_start: int
    ref_end: int
    q_start: int
    q_end: int
    raw_type: str
    subtype: str
    blocks: list = dc_field(default_factory=list, repr=False)
    true_q_len: int = 0
    sub_spans: list = dc_field(default_factory=list, repr=False)
    cluster_spans: Tuple[str, int, int, int, int, str] = ()

    @property
    def ref_len(self) -> int:
        return self.ref_end - self.ref_start

    @property
    def q_len(self) -> int:
        """Span-based query length (max qe - min qs). May include inter-block gaps."""
        return self.q_end - self.q_start

    @property
    def is_inter_chr(self) -> bool:
        return self.raw_type == 'BND'

    @property
    def is_inverted(self) -> bool:
        return self.subtype.startswith('Inv')

    @property
    def is_fragmented(self) -> bool:
        """True when sub-blocks have ref gaps > 100 bp between them (non-contiguous source)."""
        if len(self.sub_spans) <= 1:
            return False
        return any((self.sub_spans[k + 1][0] - self.sub_spans[k][1] > 100 for k in range(len(self.sub_spans) - 1)))

def _classify_orphan_cluster(c_chr: str, c_strand: str, c_r_start: int, c_r_end: int, anchor_chr: str, anchor_strand: str, anchor_bp_left: int, anchor_bp_right: int, tandem_max_dist: int, proximal_max_dist: int) -> tuple:
    inv_prefix = 'Inv' if c_strand != anchor_strand else ''
    if c_chr != anchor_chr:
        return ('BND', f'{inv_prefix}DUP:DISP:INTER')
    gap_left = min(anchor_bp_left, anchor_bp_right)
    gap_right = max(anchor_bp_left, anchor_bp_right)
    dist_to_gap = max(0, gap_left - c_r_end, c_r_start - gap_right)
    orphan_len = max(1, c_r_end - c_r_start)
    is_strictly_inside = c_r_start >= gap_left and c_r_end <= gap_right
    if is_strictly_inside:
        if c_strand != anchor_strand:
            return ('CSV', 'InvDUP')
        else:
            return ('CSV', 'LOCAL_RETAIN')
    if dist_to_gap <= tandem_max_dist and dist_to_gap <= orphan_len:
        return ('CSV', f'{inv_prefix}DUP:TANDEM')
    elif dist_to_gap <= proximal_max_dist:
        return ('CSV', f'{inv_prefix}DUP:PROXIMAL')
    else:
        return ('CSV', f'{inv_prefix}DUP:DISP:INTRA')

def condense_orphan_blocks(orphans: List[AlignedSegment], anchor_chr: str, anchor_strand: str, anchor_bp_left: int, anchor_bp_right: int, max_ref_dist: int=100000, max_query_dist: int=50000) -> List[OrphanCluster]:
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
    TANDEM_MAX = getattr(SVConfig, 'TANDEM_MAX_REF_DIST', 500)
    PROXIMAL_MAX = getattr(SVConfig, 'PROXIMAL_MAX_REF_DIST', 10000)

    def _finalize(c_chr, c_strand, c_r_start, c_r_end, c_q_start, c_q_end, c_blocks):
        (raw_type, subtype) = _classify_orphan_cluster(c_chr, c_strand, c_r_start, c_r_end, anchor_chr, anchor_strand, anchor_bp_left, anchor_bp_right, TANDEM_MAX, PROXIMAL_MAX)
        true_q_len = c_q_end - c_q_start
        sub_spans = [(b.ref_start, b.ref_end, b.qs, b.qe, b.strand) for b in c_blocks]
        cluster_spans = (c_chr, c_r_start, c_r_end, c_q_start, c_q_end, c_strand)
        return OrphanCluster(ref_name=c_chr, strand=c_strand, ref_start=c_r_start, ref_end=c_r_end, q_start=c_q_start, q_end=c_q_end, raw_type=raw_type, subtype=subtype, blocks=c_blocks, true_q_len=true_q_len, sub_spans=sub_spans, cluster_spans=cluster_spans)
    clusters = []
    seg0 = orphans[0]
    (curr_chr, curr_strand) = (seg0.ref_name, seg0.strand)
    (curr_r_start, curr_r_end) = (seg0.ref_start, seg0.ref_end)
    (curr_q_start, curr_q_end) = (seg0.qs, seg0.qe)
    curr_blocks = [seg0]
    for seg in orphans[1:]:
        ref_edge_dist = min(abs(seg.ref_start - curr_r_end), abs(seg.ref_end - curr_r_start))
        mergeable = seg.ref_name == curr_chr and seg.strand == curr_strand and (ref_edge_dist <= max_ref_dist) and (seg.qs - curr_q_end <= max_query_dist)
        if mergeable:
            curr_r_start = min(curr_r_start, seg.ref_start)
            curr_r_end = max(curr_r_end, seg.ref_end)
            curr_q_end = max(curr_q_end, seg.qe)
            curr_blocks.append(seg)
        else:
            clusters.append(_finalize(curr_chr, curr_strand, curr_r_start, curr_r_end, curr_q_start, curr_q_end, curr_blocks))
            (curr_chr, curr_strand) = (seg.ref_name, seg.strand)
            (curr_r_start, curr_r_end) = (seg.ref_start, seg.ref_end)
            (curr_q_start, curr_q_end) = (seg.qs, seg.qe)
            curr_blocks = [seg]
    clusters.append(_finalize(curr_chr, curr_strand, curr_r_start, curr_r_end, curr_q_start, curr_q_end, curr_blocks))
    return clusters

def resolve_chain_svs(chains: List['AlignedSegment'], main_strand: str, main_idx_set: Set[int], read_name: str, query_seq: Optional[str]) -> List['SVCandidate']:
    """Master controller that extracts linear/complex SVs from alignment chain transitions (Strict Backbone)."""
    candidates = []
    complex_cluster_idx = 0
    clean_block_id = 0
    if not main_idx_set:
        return candidates
    main_idx = sorted(list(main_idx_set))

    def extract_insertion_seq(qs: int, qe: int, strand: str) -> str:
        if not query_seq:
            return ''
        seq_fwd = query_seq[qs:qe]
        return seq_fwd.translate(RC_TRANS)[::-1] if strand == '-' else seq_fwd

    def check_and_extract_linear_sv(b1: 'AlignedSegment', b2: 'AlignedSegment', chrom: str, strand: str, source_label: str, isbackbone=True, extra_info: dict=None) -> Optional['SVCandidate']:
        if b1.strand != b2.strand or b1.ref_name != b2.ref_name:
            return None
        dist_q = b2.qs - b1.qe
        if strand == '+':
            (ref_break_1, ref_break_2) = (b1.ref_end, b2.ref_start)
        else:
            (ref_break_1, ref_break_2) = (b1.ref_start, b2.ref_end)
        dist_r = ref_break_2 - ref_break_1 if strand == '+' else ref_break_1 - ref_break_2
        (left_bp, right_bp) = (min(ref_break_1, ref_break_2), max(ref_break_1, ref_break_2))
        sv_len = abs(dist_q - dist_r)
        base_info = {'SOURCE': source_label}
        if extra_info:
            base_info.update(extra_info)
        if abs(dist_q) < getattr(SVConfig, 'GAP_TOLERANCE', 50) and dist_r < -getattr(SVConfig, 'MIN_SV_SIZE', 30):
            if 'COMPLEX' in source_label:
                (ref_s1, ref_e1) = (min(b1.ref_start, b1.ref_end), max(b1.ref_start, b1.ref_end))
                (ref_s2, ref_e2) = (min(b2.ref_start, b2.ref_end), max(b2.ref_start, b2.ref_end))
                ref_overlap = max(0, min(ref_e1, ref_e2) - max(ref_s1, ref_s2))
                if ref_overlap == 0 or abs(dist_r) > max(3 * (ref_e2 - ref_s2), 1000000):
                    return None
            dup_len = abs(dist_r) + dist_q
            filter_tag = 'OVERSIZE' if getattr(SVConfig, 'MAX_SV_SIZE', -1) != -1 and dup_len >= getattr(SVConfig, 'MAX_SV_SIZE', -1) else None
            sv_type = 'INS' if dup_len < getattr(SVConfig, 'TANDEM_DUP_INS_THRESHOLD', 1000) else 'DUP'
            pos = right_bp if sv_type == 'INS' else left_bp
            if sv_type == 'INS':
                base_info['TANDEM_SOURCE'] = f'{chrom}:{left_bp}-{right_bp}'
                if query_seq:
                    ins_qs = min(b1.qe, b2.qs)
                    ins_qe = max(0, min(ins_qs + dup_len, len(query_seq)))
                    if ins_qe > ins_qs:
                        base_info['SEQ'] = extract_insertion_seq(ins_qs, ins_qe, strand)
            else:
                base_info['END'] = right_bp
            base_info['SUBTYPE'] = 'DUP:TANDEM'
            cand = SVCandidate(sv_type, chrom, pos, dup_len, read_name, info=base_info, filter_str=filter_tag)
            (cand.q_start, cand.q_end) = (min(b1.qe, b2.qs), min(b1.qe, b2.qs) + dup_len)
            return cand
        if sv_len < 20:
            return None
        if sv_len >= getattr(SVConfig, 'MIN_SV_SIZE', 30) + getattr(SVConfig, 'INSDEL_TOLERANCE', 10):
            if abs(dist_q) < getattr(SVConfig, 'GAP_TOLERANCE', 50) or abs(dist_r) < getattr(SVConfig, 'GAP_TOLERANCE', 50):
                if abs(dist_r) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000):
                    return None
                sv_type = 'DEL' if dist_r > dist_q else 'INS'
                if sv_type == 'INS' and (getattr(SVConfig, 'ENABLE_LARGE_SEQUENCE', False) or sv_len < getattr(SVConfig, 'MAX_ALLELE_LEN', 5000)):
                    if query_seq:
                        (ins_qs, ins_qe) = (min(b1.qe, b2.qs), max(0, min(min(b1.qe, b2.qs) + sv_len, len(query_seq))))
                        base_info['SEQ'] = extract_insertion_seq(ins_qs, ins_qe, strand)
                cand = SVCandidate(sv_type, chrom, left_bp, sv_len, read_name, info=base_info)
                (cand.q_start, cand.q_end) = (min(b1.qe, b2.qs), max(b1.qe, b2.qs))
                # Preserve the junction geometry separately from its net SVLEN.
                # A DEL net length excludes query filler; the reference gap does
                # not. Edge-CIGAR provenance needs the actual consumed intervals.
                cand.reference_gap=(left_bp,right_bp)
                return cand
            else:
                if abs(dist_r) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000):
                    return None
                if isbackbone and dist_q > 0 and (dist_r > 0):
                    (ins_len, del_len) = (int(dist_q), int(dist_r))
                    if ins_len >= getattr(SVConfig, 'GAP_TOLERANCE', 50) and del_len >= getattr(SVConfig, 'GAP_TOLERANCE', 50):
                        if sv_len >= getattr(SVConfig, 'MIN_CSV_SIZE', getattr(SVConfig, 'MIN_SV_SIZE', 30) * 1.2):
                            base_info.update({'SUBTYPE': 'INS_DEL', 'INS_LEN': ins_len, 'DEL_LEN': del_len, 'END': left_bp + del_len})
                            if query_seq:
                                (ins_qs, ins_qe) = (min(b1.qe, b2.qs), max(0, min(min(b1.qe, b2.qs) + ins_len, len(query_seq))))
                                base_info['SEQ'] = extract_insertion_seq(ins_qs, ins_qe, strand)
                            cand = SVCandidate('CSV', chrom, left_bp, sv_len, read_name, info=base_info)
                            (cand.q_start, cand.q_end) = (min(b1.qe, b2.qs), max(b1.qe, b2.qs))
                            return cand
        return None

    def _extract_and_merge_nested_nodes(nodes: List['AlignedSegment'], source_label: str, macro_event_id: str, nested_dict: dict):
        """提取内部嵌套微变异，统一绑定宏观 EVENT ID，并动态分配唯一 Candidate ID"""
        for k in range(len(nodes) - 1):
            (n1, n2) = (nodes[k], nodes[k + 1])
            if n1.ref_name != n2.ref_name:
                continue
            internal_cand = check_and_extract_linear_sv(b1=n1, b2=n2, chrom=n1.ref_name, strand=n1.strand, source_label=source_label, isbackbone=False, extra_info={'EVENT': macro_event_id})
            if not internal_cand:
                continue
            merge_key = next((_mk for _mk in nested_dict if internal_cand.raw_type == _mk[0] and abs(internal_cand.pos - _mk[1]) <= 10), None)
            if merge_key is None:
                nest_idx = len(nested_dict)
                internal_cand.id = f'{macro_event_id}_NEST_{nest_idx}'
                internal_cand.info['MERGED_BLOCKS'] = 1
                nested_dict[internal_cand.raw_type, internal_cand.pos] = internal_cand
            else:
                existing = nested_dict[merge_key]
                e_span = existing.q_end - existing.q_start
                n_span = internal_cand.q_end - internal_cand.q_start
                is_same_detection = True
                if e_span > 0 and n_span > 0:
                    overlap = max(0, min(existing.q_end, internal_cand.q_end) - max(existing.q_start, internal_cand.q_start))
                    shorter = max(1, min(e_span, n_span))
                    is_same_detection = overlap / shorter > 0.5
                if is_same_detection:
                    if internal_cand.sv_len > existing.sv_len:
                        existing.sv_len = internal_cand.sv_len
                        for key in ['INS_LEN', 'DEL_LEN', 'SEQ', 'END']:
                            if key in internal_cand.info:
                                existing.info[key] = internal_cand.info[key]
                        existing.end = internal_cand.end
                else:
                    existing.sv_len += internal_cand.sv_len
                    existing.info['MERGED_BLOCKS'] = existing.info.get('MERGED_BLOCKS', 1) + 1
                    existing.info['INS_LEN'] = existing.info.get('INS_LEN', 0) + internal_cand.info.get('INS_LEN', 0)
                    existing.info['DEL_LEN'] = existing.info.get('DEL_LEN', 0) + internal_cand.info.get('DEL_LEN', 0)
                    if internal_cand.info.get('SUBTYPE') == 'DUP:TANDEM':
                        existing.info['SUBTYPE'] = 'MULTI_DUP:TANDEM'
                    if 'SEQ' in existing.info and 'SEQ' in internal_cand.info:
                        existing.info['SEQ'] += internal_cand.info['SEQ']
                    if 'END' in internal_cand.info and 'END' in existing.info:
                        existing.info['END'] = max(existing.info['END'], internal_cand.info['END'])
                    existing.q_end = max(existing.q_end, internal_cand.q_end)
                    if 'END' in existing.info:
                        existing.end = int(existing.info['END'])
                    elif existing.type in ('BND', 'TRA'):
                        existing.end = existing.pos
                    elif existing.type == 'INS':
                        existing.end = existing.pos + int(existing.info.get('REPLACED_REF_LEN', 0))
                    else:
                        existing.end = existing.pos + existing.sv_len

    def process_edge_orphans(bb_node: 'AlignedSegment', orphans: List['AlignedSegment'], is_head: bool):
        if not orphans:
            return
        bb_pos = bb_node.ref_start if bb_node.strand == '+' else bb_node.ref_end
        if not is_head:
            bb_pos = bb_node.ref_end if bb_node.strand == '+' else bb_node.ref_start
        (anchor_chr, anchor_strand) = (bb_node.ref_name, bb_node.strand)
        global_event_id = f"{read_name}_{('HEAD' if is_head else 'TAIL')}_EDGE"
        merged_clusters = condense_orphan_blocks(orphans, anchor_chr=anchor_chr, anchor_strand=anchor_strand, anchor_bp_left=bb_pos, anchor_bp_right=bb_pos)
        nested_merged_dict = {}
        _extract_and_merge_nested_nodes(orphans, 'COMPLEX_EDGE', global_event_id, nested_merged_dict)
        candidates.extend(nested_merged_dict.values())
        total_ins_len = max((o.qe for o in orphans)) - min((o.qs for o in orphans))
        sub_type_parts = []
        for (idx, m) in enumerate(merged_clusters):
            if len(merged_clusters) == 1 and m.strand != anchor_strand and (m.ref_name == anchor_chr):
                inv_span_len = m.ref_end - m.ref_start
                TRA_GAP_MIN = getattr(SVConfig, 'MIN_TRA_SIZE', 50000)
                TRA_RATIO_THRESHOLD = getattr(SVConfig, 'TRA_RATIO_THRESHOLD', 5.0)
                raw_gap = abs(m.ref_start - bb_pos) if is_head else abs(bb_pos - m.ref_end)
                if raw_gap > TRA_GAP_MIN or raw_gap > inv_span_len * TRA_RATIO_THRESHOLD:
                    m.subtype = m.subtype.replace('TANDEM', 'DISP:INTRA')
                else:
                    m.subtype = 'INV'
            sub_type_parts.append(m.subtype)
            span_info = {'EVENT': global_event_id, 'SRCREG': f'{m.ref_name}:{m.ref_start}-{m.ref_end}', 'SOURCE': 'COMPLEX_EDGE', 'SUBTYPE': m.subtype, 'DEL_LEN': 0, 'INS_LEN': m.true_q_len, 'END': bb_pos}
            if m.is_fragmented:
                span_info['SUBREGS'] = ','.join((f'{s}-{e}' for (s, e, qs, qe, strd) in m.sub_spans))
            span_cand = SVCandidate('CSV', anchor_chr, bb_pos, m.true_q_len, read_name, info=span_info, filter_str='EDGE_COMPLEX')
            span_cand.id = f'{global_event_id}_SPAN_{idx}'
            candidates.append(span_cand)
        eventtype = merged_clusters[0].subtype if len(merged_clusters) == 1 else 'COMPLEX_' + '_'.join(sub_type_parts)
        summary_cand = SVCandidate('CSV', anchor_chr, bb_pos, total_ins_len, read_name, info={'EVENT': global_event_id, 'EVENTTYPE': eventtype, 'SUBTYPE': 'COMPLEX', 'DEL_LEN': 0, 'INS_LEN': total_ins_len, 'ORPHAN_CLUSTERS': len(merged_clusters)}, filter_str='EDGE_COMPLEX')
        summary_cand.id = f'{global_event_id}_SUMMARY'
        candidates.append(summary_cand)
        closest_orphan = orphans[-1] if is_head else orphans[0]
        orphan_pos = closest_orphan.ref_end if closest_orphan.strand == '+' else closest_orphan.ref_start
        if not is_head:
            orphan_pos = closest_orphan.ref_start if closest_orphan.strand == '+' else closest_orphan.ref_end
        strands = closest_orphan.strand + anchor_strand if is_head else anchor_strand + closest_orphan.strand
        tra_type = 'INTER_CHR' if anchor_chr != closest_orphan.ref_name else 'INTRA_CHR'
        c1 = SVCandidate('BND', anchor_chr, bb_pos, 0, read_name, info={'EVENT': global_event_id, 'CHR2': closest_orphan.ref_name, 'POS2': orphan_pos, 'TRA_TYPE': tra_type, 'STRANDS': strands, 'SOURCE': 'COMPLEX_EDGE'}, filter_str='EDGE_COMPLEX')
        c2 = SVCandidate('BND', closest_orphan.ref_name, orphan_pos, 0, read_name, info={'EVENT': global_event_id, 'CHR2': anchor_chr, 'POS2': bb_pos, 'TRA_TYPE': tra_type, 'STRANDS': strands[::-1], 'SOURCE': 'COMPLEX_EDGE'}, filter_str='EDGE_COMPLEX')
        c1.id = f'{global_event_id}_BND_EDGE_B'
        c2.id = f'{global_event_id}_BND_EDGE_A'
        c1.mateid = c2.id
        c2.mateid = c1.id
        c1.info['MATEID'] = c1.mateid
        c2.info['MATEID'] = c2.mateid
        candidates.extend([c1, c2])
    if main_idx[0] > 0:
        process_edge_orphans(bb_node=chains[main_idx[0]], orphans=chains[:main_idx[0]], is_head=True)
    if main_idx[-1] < len(chains) - 1:
        process_edge_orphans(bb_node=chains[main_idx[-1]], orphans=chains[main_idx[-1] + 1:], is_head=False)
    for i in range(len(main_idx) - 1):
        (idx_L, idx_R) = (main_idx[i], main_idx[i + 1])
        (bb_left, bb_right) = (chains[idx_L], chains[idx_R])
        if not (bb_left.is_backbone and bb_right.is_backbone):
            continue
        orphans_in_between = chains[idx_L + 1:idx_R]
        (anchor_chr, anchor_strand) = (bb_left.ref_name, bb_left.strand)
        anchor_bp_left = bb_left.ref_end if anchor_strand == '+' else bb_left.ref_start
        anchor_bp_right = bb_right.ref_start if anchor_strand == '+' else bb_right.ref_end
        (anchor_pos, anchor_end) = (anchor_bp_left, anchor_bp_right)
        dist_q = bb_right.qs - bb_left.qe
        dist_r = bb_right.ref_start - bb_left.ref_end if anchor_strand == '+' else bb_left.ref_start - bb_right.ref_end
        if orphans_in_between:
            clean_block_id += 1
            merged_clusters = condense_orphan_blocks(orphans_in_between, anchor_chr=anchor_chr, anchor_strand=anchor_strand, anchor_bp_left=anchor_bp_left, anchor_bp_right=anchor_bp_right)
            complex_cluster_idx += 1
            global_event_id = f'{read_name}_CPLX_{complex_cluster_idx}'
            if len(merged_clusters) == 1:
                m = merged_clusters[0]
                is_inv = m.strand != anchor_strand
                is_same_chr = m.ref_name == anchor_chr
                if is_inv and is_same_chr:
                    (inv_start, inv_end) = (m.ref_start, m.ref_end)
                    inv_span_len = inv_end - inv_start
                    if anchor_strand == '+':
                        (pos_outer_left, pos_outer_right) = (bb_left.ref_end, bb_right.ref_start)
                        (q_raw_left, q_raw_right) = (m.q_start - bb_left.qe, bb_right.qs - m.q_end)
                    else:
                        (pos_outer_left, pos_outer_right) = (bb_right.ref_end, bb_left.ref_start)
                        (q_raw_left, q_raw_right) = (bb_right.qs - m.q_end, m.q_start - bb_left.qe)
                    (raw_gap_left, raw_gap_right) = (inv_start - pos_outer_left, pos_outer_right - inv_end)
                    (delta_left, delta_right) = (q_raw_left - raw_gap_left, q_raw_right - raw_gap_right)
                    (true_ins_left, true_del_left) = (max(0, delta_left), max(0, -delta_left))
                    (true_ins_right, true_del_right) = (max(0, delta_right), max(0, -delta_right))
                    total_span = pos_outer_right - pos_outer_left
                    TRA_GAP_MIN = getattr(SVConfig, 'MIN_TRA_SIZE', 50000)
                    TRA_RATIO_THRESHOLD = getattr(SVConfig, 'TRA_RATIO_THRESHOLD', 5.0)
                    is_remote_inv = abs(raw_gap_left) > TRA_GAP_MIN or abs(raw_gap_right) > TRA_GAP_MIN or abs(raw_gap_left) > inv_span_len * TRA_RATIO_THRESHOLD or (abs(raw_gap_right) > inv_span_len * TRA_RATIO_THRESHOLD)
                    if is_remote_inv:
                        m.subtype = m.subtype.replace('TANDEM', 'DISP:INTRA')
                    else:
                        MIN_SV = getattr(SVConfig, 'MIN_SV_SIZE', 30)
                        total_ins = true_ins_left + true_ins_right
                        total_del = true_del_left + true_del_right
                        is_ins = total_ins > MIN_SV
                        is_del = total_del > MIN_SV
                        if is_ins and is_del:
                            inv_subtype = 'INV-INS-DEL'
                        elif is_ins:
                            inv_subtype = 'INV-INS'
                        elif is_del:
                            inv_subtype = 'INV-DEL'
                        else:
                            inv_subtype = 'INV'
                        info = {'EVENT': global_event_id, 'SUBTYPE': inv_subtype, 'INV_SPAN': int(inv_span_len), 'BP_DETAILS': f'LC:{int(q_raw_left)},RC:{int(q_raw_right)},LF:{int(raw_gap_left)},RF:{int(raw_gap_right)}'}
                        if is_ins:
                            info['BP_INS'] = int(delta_left + delta_right)
                        if is_del:
                            info['BP_DEL'] = int(-(delta_right + delta_left))
                        if inv_subtype == 'INV':
                            inv_cand = SVCandidate('INV', anchor_chr, inv_start, inv_span_len, read_name, info=info)
                            inv_cand.id = f'{global_event_id}_INV'
                            candidates.append(inv_cand)
                        else:
                            macro_ins_len = int(q_raw_left + inv_span_len + q_raw_right)
                            info.update({'DEL_LEN': int(total_span), 'INS_LEN': macro_ins_len, 'END': int(pos_outer_right)})
                            csv_cand = SVCandidate('CSV', anchor_chr, pos_outer_left, macro_ins_len, read_name, info=info)
                            csv_cand.id = f'{global_event_id}_CSV'
                            candidates.append(csv_cand)
                        if len(m.blocks) > 1:
                            nested_merged_dict = {}
                            _extract_and_merge_nested_nodes(m.blocks, 'COMPLEX', global_event_id, nested_merged_dict)
                            candidates.extend(nested_merged_dict.values())
                        continue
            nested_merged_dict = {}
            all_nodes = [bb_left] + orphans_in_between + [bb_right]
            _extract_and_merge_nested_nodes(all_nodes, 'COMPLEX', global_event_id, nested_merged_dict)
            candidates.extend(nested_merged_dict.values())
            sub_type_info = []
            sub_type_info.append(f'backbone,{bb_left.ref_name},{anchor_bp_left},{anchor_bp_right},{bb_left.qe},{bb_right.qs},{main_strand}')
            sub_type_parts = []
            (prev_chr, prev_strand, prev_pos) = (anchor_chr, anchor_strand, anchor_pos)
            for (idx, m) in enumerate(merged_clusters):
                sub_type_parts.append(m.subtype)
                sub_type_info.append(f'{m.subtype},{m.cluster_spans[0]},{m.cluster_spans[1]},{m.cluster_spans[2]},{m.cluster_spans[3]},{m.cluster_spans[4]},{m.cluster_spans[5]}')
                m_enter_pos = m.ref_start if m.strand == '+' else m.ref_end
                m_leave_pos = m.ref_end if m.strand == '+' else m.ref_start
                is_inter_chr = prev_chr != m.ref_name
                is_distant_intra = not is_inter_chr and abs(prev_pos - m_enter_pos) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000)
                is_distant_jump = is_inter_chr or is_distant_intra
                bnd_entry = SVCandidate('BND', prev_chr, prev_pos, 0, read_name, info={'EVENT': global_event_id, 'CHR2': m.ref_name, 'POS2': m_enter_pos, 'TRA_TYPE': 'INTER_CHR' if is_inter_chr else 'INTRA_CHR', 'STRANDS': prev_strand + m.strand, 'SOURCE': 'COMPLEX'}, filter_str=None if is_distant_jump else 'LOCAL_BND')
                bnd_entry.id = f'{global_event_id}_BND_{idx}_A'
                bnd_entry.mateid = f'{global_event_id}_BND_{idx}_B'
                bnd_exit = SVCandidate('BND', m.ref_name, m_enter_pos, 0, read_name, info={'EVENT': global_event_id, 'CHR2': prev_chr, 'POS2': prev_pos, 'TRA_TYPE': 'INTER_CHR' if is_inter_chr else 'INTRA_CHR', 'STRANDS': m.strand + prev_strand, 'SOURCE': 'COMPLEX'}, filter_str=None if is_distant_jump else 'LOCAL_BND')
                bnd_exit.id = bnd_entry.mateid
                bnd_exit.mateid = bnd_entry.id
                bnd_entry.info['MATEID'] = bnd_entry.mateid
                bnd_exit.info['MATEID'] = bnd_exit.mateid
                candidates.extend([bnd_entry, bnd_exit])
                (prev_chr, prev_strand, prev_pos) = (m.ref_name, m.strand, m_leave_pos)
                span_info = {'EVENT': global_event_id, 'SRCREG': f'{m.ref_name}:{m.ref_start}-{m.ref_end}', 'SUBTYPE': m.subtype, 'DEL_LEN': 0, 'INS_LEN': m.true_q_len}
                if m.is_fragmented:
                    span_info['SUBREGS'] = ','.join((f'{s}-{e}' for (s, e, qs, qe, strd) in m.sub_spans))
                if m.subtype != 'DUP:TANDEM':
                    if m.subtype == 'LOCAL_RETAIN':
                        pass
                    elif m.subtype == 'InvDUP':
                        span_info['END'] = m.ref_end
                        inv_len = m.ref_end - m.ref_start
                        span_cand = SVCandidate('CSV', m.ref_name, m.ref_start, inv_len, read_name, info=span_info)
                        span_cand.id = f'{global_event_id}_SPAN_{idx}'
                        candidates.append(span_cand)
                    else:
                        span_info['END'] = anchor_bp_left
                        span_cand = SVCandidate('CSV', anchor_chr, anchor_bp_left, m.true_q_len, read_name, info=span_info)
                        span_cand.id = f'{global_event_id}_SPAN_{idx}'
                        candidates.append(span_cand)
            bb_right_enter = bb_right.ref_start if bb_right.strand == '+' else bb_right.ref_end
            tail_inter = prev_chr != bb_right.ref_name
            tail_distant = not tail_inter and abs(prev_pos - bb_right_enter) > getattr(SVConfig, 'MIN_TRA_SIZE', 50000)
            if tail_inter or tail_distant:
                t1 = SVCandidate('BND', prev_chr, prev_pos, 0, read_name, info={'EVENT': global_event_id, 'CHR2': bb_right.ref_name, 'POS2': bb_right_enter, 'TRA_TYPE': 'INTER_CHR' if tail_inter else 'INTRA_CHR', 'STRANDS': prev_strand + bb_right.strand})
                t1.id = f'{global_event_id}_BND_TAIL_A'
                t1.mateid = f'{global_event_id}_BND_TAIL_B'
                t2 = SVCandidate('BND', bb_right.ref_name, bb_right_enter, 0, read_name, info={'EVENT': global_event_id, 'CHR2': prev_chr, 'POS2': prev_pos, 'TRA_TYPE': 'INTER_CHR' if tail_inter else 'INTRA_CHR', 'STRANDS': bb_right.strand + prev_strand})
                t2.id = t1.mateid
                t2.mateid = t1.id
                t1.info['MATEID'] = t1.mateid
                t2.info['MATEID'] = t2.mateid
                candidates.extend([t1, t2])
            if len(merged_clusters) >= 2:
                up_down = bool(anchor_pos > anchor_end)
                affected_ref = abs(dist_r)
                eventtype = merged_clusters[0].subtype if len(merged_clusters) == 1 else 'COMPLEX_' + '_'.join(sub_type_parts)
                summary_pos = anchor_end if up_down else anchor_pos
                summary_end = summary_pos + affected_ref
                serialized_sub_svs = '|'.join(sub_type_info)
                summary_info = {'EVENT': global_event_id, 'EVENTTYPE': eventtype, 'SUBTYPE': 'COMPLEX', 'DEL_LEN': int(affected_ref), 'INS_LEN': int(max(0, dist_q)), 'END': int(summary_end), 'ORPHAN_CLUSTERS': len(merged_clusters), 'SUB_SVS': serialized_sub_svs}
                summary_cand = SVCandidate('CSV', anchor_chr, summary_pos, max(0, dist_q), read_name, info=summary_info, filter_str='COMPLEX')
                summary_cand.id = f'{global_event_id}_SUMMARY'
                candidates.append(summary_cand)
        else:
            backbone_cand = check_and_extract_linear_sv(b1=bb_left, b2=bb_right, chrom=anchor_chr, strand=anchor_strand, source_label='backbone')
            if backbone_cand:
                candidates.append(backbone_cand)
    return candidates
