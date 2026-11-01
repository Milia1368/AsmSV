"""Original SSV gap kernel, lifted verbatim from the legacy detector.

The frozen closure retains legacy behavior for comparisons. The structured
detector intercepts duplicated reference visits before this closure, so its old
tandem branch is not the production duplication rule. Non-tandem INS/DEL and
sequence slicing remain unchanged. Sequence DELINS is handled outside this kernel.
"""
from typing import Optional
from asmsv_next.alignment import AlignedSegment
from asmsv_next.config import SVConfig
from asmsv_next.structures import SVCandidate
RC_TRANS = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')


def make_gap_caller(read_name, query_seq):
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
    return check_and_extract_linear_sv
