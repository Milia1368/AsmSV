"""Frozen 0.5.0 parser oracle; never imported by the caller."""
from typing import Optional, List, Dict
from asmsv_next.alignment import AlignedSegment, RC_TRANS

class LegacySegment(AlignedSegment):
    def parse_intra_alignment_sv(self, query_seq: Optional[str]=None, min_size: int=50, merge_gap_tol: int=50, *, original=False, _geometry=None) -> List[Dict]:
        """Parses INDELs found directly inside the CIGAR string."""
        svs = []
        raw=self.score_root.raw
        cigar=raw.cigar if original else self.cigar
        ref_start=raw.ref_start if original else self.ref_start
        qs,qe=(raw.qs,raw.qe) if original else (self.qs,self.qe)
        if _geometry is not None:cigar,ref_start,qs,qe=_geometry
        if not cigar:
            return svs
        curr_ref, curr_bam_seq_idx = (ref_start, 0)
        left_clip = 0
        for op, length in cigar:
            if op in {4, 5}:
                left_clip += length
            elif op in {0, 1, 7, 8}:
                break
        p_ins, ins_gap = (None, 0)
        p_del, del_gap = (None, 0)

        def flush_ins():
            nonlocal p_ins
            if p_ins and p_ins['len'] >= min_size:
                svs.append({'type': 'INS', 'chr': self.ref_name, 'pos': p_ins['pos'], 'len': p_ins['len'], 'read': self.q_name, 'info': {'SEQ': p_ins['seq'], 'SOURCE': 'CIGAR_MERGED' if 'MERGED' in p_ins.get('flag', '') else 'CIGAR'}})
            p_ins = None

        def flush_del():
            nonlocal p_del
            if p_del and p_del['len'] >= min_size:
                svs.append({'type': 'DEL', 'chr': self.ref_name, 'pos': p_del['pos'], 'len': p_del['len'], 'read': self.q_name, 'info': {'SOURCE': 'CIGAR_MERGED' if 'MERGED' in p_del.get('flag', '') else 'CIGAR'}})
            p_del = None
        for op, length in cigar:
            consumes_seq = op in {0, 1, 4, 5, 7, 8}
            consumes_ref = op in {0, 2, 3, 7, 8}
            if op == 1:
                ins_seq = ''
                if query_seq:
                    if not self.is_reverse:
                        g_start = qs + (curr_bam_seq_idx - left_clip)
                        g_end = g_start + length
                        if 0 <= g_start < g_end <= len(query_seq):
                            ins_seq = query_seq[g_start:g_end]
                    else:
                        g_end = qe - (curr_bam_seq_idx - left_clip)
                        g_start = g_end - length
                        if 0 <= g_start < g_end <= len(query_seq):
                            ins_seq = query_seq[g_start:g_end].translate(RC_TRANS)[::-1]
                sv_pos = curr_ref
                if p_ins:
                    if ins_gap <= merge_gap_tol:
                        p_ins['len'] += length
                        p_ins['seq'] += ins_seq
                        p_ins['flag'] = 'MERGED'
                        ins_gap = 0
                    else:
                        flush_ins()
                        p_ins = {'pos': sv_pos, 'len': length, 'seq': ins_seq}
                        ins_gap = 0
                else:
                    p_ins = {'pos': sv_pos, 'len': length, 'seq': ins_seq}
                    ins_gap = 0
                if p_del:
                    del_gap += length
            elif op == 2:
                sv_pos = curr_ref
                if p_del:
                    if del_gap <= merge_gap_tol:
                        p_del['len'] += length
                        p_del['flag'] = 'MERGED'
                        del_gap = 0
                    else:
                        flush_del()
                        p_del = {'pos': sv_pos, 'len': length}
                        del_gap = 0
                else:
                    p_del = {'pos': sv_pos, 'len': length}
                    del_gap = 0
                if p_ins:
                    ins_gap += length
            else:
                if consumes_ref and p_ins:
                    ins_gap += length
                if consumes_seq and consumes_ref and p_del:
                    del_gap += length
            if p_ins and ins_gap > merge_gap_tol:
                flush_ins()
            if p_del and del_gap > merge_gap_tol:
                flush_del()
            if consumes_seq:
                curr_bam_seq_idx += length
            if consumes_ref:
                curr_ref += length
        flush_ins()
        flush_del()
        return svs
