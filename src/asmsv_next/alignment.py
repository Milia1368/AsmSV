"""One original alignment block, immutable evidence root and idempotent views."""
from asmsv_next.aq import RawAlignment, EvidenceRoot
from asmsv_next.config import SVConfig
from asmsv_next.constants import REWARD_FLOOR
from typing import Optional, List, Dict

RC_TRANS=str.maketrans('ATCGNatcgn','TAGCNtagcn')


class AlignedSegment:
    __slots__=('q_name','flag','cigar','ref_name','ref_start','ref_end','is_reverse','mapq',
        'as_score','strand','qs','qe','density','is_valid','is_backbone','score_root',
        'summary','aq','raw_as','aq_sequence','mapq_factor','rejection_reasons','quality_flags','reward_rescued')

    def __init__(self,from_dict):
        data=from_dict
        self.q_name=data['q_name'];self.flag=data.get('flag',0)
        self.cigar=[tuple(x) for x in data['cigar']]
        self.ref_name=data['ref_name'];self.ref_start=data['ref_start'];self.ref_end=data['ref_end']
        self.is_reverse=data.get('is_reverse',False);self.mapq=data.get('mapq',0)
        self.strand=data.get('strand','-' if self.is_reverse else '+')
        self.score_root=data.get('_score_root')
        if self.score_root is None:
            self.score_root=EvidenceRoot(RawAlignment.from_dict(data),SVConfig.CALLER_PRESET.aq)
        self.qs=data.get('qs',self.score_root.raw.qs);self.qe=data.get('qe',self.score_root.raw.qe)
        self.raw_as=self.score_root.raw.raw_as;self.is_backbone=False
        self._apply_score()
        if '_score_root' not in data:self._refine_edges()

    def _apply_score(self):
        p=SVConfig.CALLER_PRESET;self.summary=self.score_root.index.score(self.qs,self.qe)
        self.aq_sequence=self.summary.aq;self.mapq_factor=p.aq.mapq_factor(self.mapq)
        self.aq=self.aq_sequence*self.mapq_factor if self.aq_sequence is not None and self.mapq_factor is not None else None
        reward=self.summary.node_reward
        self.reward_rescued=reward is not None and reward<REWARD_FLOOR
        reward=max(REWARD_FLOOR,reward) if reward is not None else None
        self.as_score=reward*self.mapq_factor if reward is not None and self.mapq_factor is not None else 0.0
        self.density=self.as_score/(self.qe-self.qs) if self.qe>self.qs else 0.0
        # Native AS is diagnostic only, including malformed/nonfinite values.
        reasons=[x for x in self.score_root.issues if x!='INVALID_NATIVE_AS']
        if not SVConfig.MIN_MAPQ<=self.mapq<255:reasons.append('MAPQ_INELIGIBLE')
        if self.aq is None:reasons.append('UNKNOWN_DIFFERENCES')
        if self.as_score<=0:reasons.append('NONPOSITIVE_REWARD')
        if self.qe<=self.qs:reasons.append('EMPTY_QUERY_INTERVAL')
        if not self.summary.positional and (self.qs,self.qe)!=(self.score_root.raw.qs,self.score_root.raw.qe):
            reasons.append('UNLOCATED_SLICE_DIFFERENCES')
        self.rejection_reasons=tuple(sorted(set(reasons)));self.is_valid=not reasons
        self.quality_flags=tuple(x for x in (
            'INVALID_NATIVE_AS' if 'INVALID_NATIVE_AS' in self.score_root.issues else None,
            'REWARD_FLOORED' if self.reward_rescued else None) if x)

    def _refine_edges(self):
        if not self.summary.positional:return
        anchors=[r for r in self.score_root.index.local_runs(self.qs,self.qe)
                 if r.op==7 and min(self.qe,r.qe)-max(self.qs,r.qs)>=SVConfig.CALLER_PRESET.aq.min_anchor]
        if not anchors:return
        start,end=max(self.qs,anchors[0].qs),min(self.qe,anchors[-1].qe)
        if (start,end)==(self.qs,self.qe):return
        cigar,rs,re=self.score_root.cigar_view(start,end)
        self.qs,self.qe,self.cigar,self.ref_start,self.ref_end=start,end,cigar,rs,re
        self._apply_score()

    def split_at_query_coords(self,start,end):
        start,end=max(start,self.qs),min(end,self.qe)
        bounds=self.score_root.mapped_bounds(start,end)
        if bounds is None:return None
        start,end=bounds
        cigar,rs,re=self.score_root.cigar_view(start,end)
        return AlignedSegment(dict(q_name=self.q_name,flag=self.flag,cigar=cigar,
            ref_name=self.ref_name,ref_start=rs,ref_end=re,is_reverse=self.is_reverse,
            mapq=self.mapq,strand=self.strand,qs=start,qe=end,_score_root=self.score_root))

    def to_block_tuple(self):
        return self.qs,self.qe,self.ref_start,self.ref_end,self.ref_name,self.strand

    def parse_intra_alignment_sv(self, query_seq: Optional[str]=None, min_size: int=50, merge_gap_tol: int=50, *, original=False, _geometry=None, with_provenance=False) -> List[Dict]:
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
        provenance_query_idx = 0
        left_clip = 0
        for op, length in cigar:
            if op in {4, 5}:
                left_clip += length
            elif op in {0, 1, 7, 8}:
                break
        p_ins, ins_gap = (None, 0)
        p_del, del_gap = (None, 0)

        def attach_provenance(pending):
            if with_provenance:
                svs[-1]['evidence'] = dict(alignment_id=raw.alignment_id,
                    view_query=[qs,qe],parts=pending.get('_parts', []))

        def flush_ins():
            nonlocal p_ins
            if p_ins and p_ins['len'] >= min_size:
                svs.append({'type': 'INS', 'chr': self.ref_name, 'pos': p_ins['pos'], 'len': p_ins['len'], 'read': self.q_name, 'info': {'SEQ': p_ins['seq'], 'SOURCE': 'CIGAR_MERGED' if 'MERGED' in p_ins.get('flag', '') else 'CIGAR'}})
                attach_provenance(p_ins)
            p_ins = None

        def flush_del():
            nonlocal p_del
            if p_del and p_del['len'] >= min_size:
                svs.append({'type': 'DEL', 'chr': self.ref_name, 'pos': p_del['pos'], 'len': p_del['len'], 'read': self.q_name, 'info': {'SOURCE': 'CIGAR_MERGED' if 'MERGED' in p_del.get('flag', '') else 'CIGAR'}})
                attach_provenance(p_del)
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
                if with_provenance:
                    point = qs + provenance_query_idx-left_clip if not self.is_reverse else qe-provenance_query_idx+left_clip
                    qstart,qend=(point,point+length) if not self.is_reverse else (point-length,point)
                    p_ins.setdefault('_parts',[]).append([1,qstart,qend,curr_ref,curr_ref])
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
                if with_provenance:
                    point = qs + provenance_query_idx-left_clip if not self.is_reverse else qe-provenance_query_idx+left_clip
                    p_del.setdefault('_parts',[]).append([2,point,point,curr_ref,curr_ref+length])
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
            if with_provenance and (consumes_seq or op == 5):
                provenance_query_idx += length
            if consumes_ref:
                curr_ref += length
        flush_ins()
        flush_del()
        return svs
