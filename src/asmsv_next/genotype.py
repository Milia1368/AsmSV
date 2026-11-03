"""Deterministic allele clustering and assembly-support genotyping."""
import logging
import json
import edlib
from typing import List, Callable, Dict, Set, Optional
from collections import defaultdict
from asmsv_next.structures import SVCandidate
from asmsv_next.config import SVConfig
logger = logging.getLogger(__name__)

def candidate_key(sv):
    """Total order for equal-position candidates, independent of worker completion."""
    metadata = {k: v for (k, v) in sv.info.items() if k not in ('SEQ', 'ALLELE_ALT_SEQ', 'REGION_ID')}
    extra = json.dumps(metadata, sort_keys=True, separators=(',', ':'), ensure_ascii=False, default=lambda x: sorted(x) if isinstance(x, set) else str(x))
    return (sv.chrom, sv.pos, sv.sv_len, sv.end, sv.type, sv.raw_type, str(sv.info.get('CHR2', '')), int(sv.info.get('POS2', 0)), str(sv.info.get('STRANDS', '')), str(sv.info.get('SEQ', '')), str(sv.info.get('ALLELE_ALT_SEQ', '')), sv.filter, sv.id, sv.mateid, tuple(sorted(sv.read_names)), tuple(sorted(sv.sources)), extra)

def calculate_seq_similarity(seq1: str, seq2: str) -> float:
    """
    Calculates sequence similarity tailored for assembly data.
    Retains physical edge trimming to mitigate aligner wobble but tightens internal edit tolerance.
    """
    if not seq1 or not seq2:
        return 0.0
    (len1, len2) = (len(seq1), len(seq2))
    trim_len = 2
    if len1 > 3 * trim_len and len2 > 3 * trim_len:
        seq1_core = seq1[trim_len:-trim_len]
        seq2_core = seq2[trim_len:-trim_len]
    else:
        (seq1_core, seq2_core) = (seq1, seq2)
    core_max_len = max(len(seq1_core), len(seq2_core))
    align_res = edlib.align(seq1_core, seq2_core, mode='NW', task='distance')
    if align_res['editDistance'] == -1:
        return 0.0
    ed = align_res['editDistance']
    BASE_EDIT_TOLERANCE = 0
    if ed <= BASE_EDIT_TOLERANCE:
        return 1.0
    adjusted_ed = ed - BASE_EDIT_TOLERANCE
    return 1.0 - adjusted_ed / core_max_len

def _match_span_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic for DEL, DUP, and INV.
    Prioritizes positional and length proximity over sequence alignment for speed.
    """
    if s1.type=='DUP' and (s1.info.get('TARGET_CHR') or s2.info.get('TARGET_CHR')):
        if s1.info.get('TARGET_CHR')!=s2.info.get('TARGET_CHR'):return False
        if abs(s1.info.get('TARGET_POS',-1)-s2.info.get('TARGET_POS',-1))>SVConfig.CLUSTER_POS_TOL:return False
        if any(s1.info.get(key)!=s2.info.get(key) for key in ('DUP_ORIENT','DUP_MODE')):return False
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
    seq1 = s1.info.get('SEQ', '')
    seq2 = s2.info.get('SEQ', '')
    if seq1 and seq2:
        sim = calculate_seq_similarity(seq1, seq2)
        return sim >= SVConfig.MIN_SEQ_SIMILARITY
    return True

def _match_bnd_tra_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic for BND and TRA (Translocations).
    Requires both breakpoints to be proximally matched.
    """
    if s1.info.get('SOURCE')=='STRUCTURE_VIEW' or s2.info.get('SOURCE')=='STRUCTURE_VIEW':
        return s1.info.get('SOURCE')==s2.info.get('SOURCE') and s1.id==s2.id
    if s1.info.get('STRANDS') != s2.info.get('STRANDS'):
        return False
    if s1.chrom != s2.chrom:
        return False
    if abs(s1.pos - s2.pos) > SVConfig.CLUSTER_POS_TOL:
        return False
    chr2_1 = s1.info.get('CHR2')
    chr2_2 = s2.info.get('CHR2')
    if chr2_1 != chr2_2:
        return False
    pos2_1 = s1.info.get('POS2', -1)
    pos2_2 = s2.info.get('POS2', -1)
    if abs(pos2_1 - pos2_2) > SVConfig.CLUSTER_POS_TOL:
        return False
    return True

def _match_composite_sv(s1: SVCandidate, s2: SVCandidate) -> bool:
    """
    Matching logic for DELINS and CPX components.
    Checks positional boundaries (POS and END) and compares INS_LEN and DEL_LEN components separately.
    """
    if s1.type=='DELINS' and s1.info.get('SEQ')!=s2.info.get('SEQ'):return False
    if s1.type=='CPX':
        if any(s1.info.get(k)!=s2.info.get(k) for k in ('CPX_TYPE','CPX_COMPONENTS','STRUCTURE_KEY')):return False
    if s1.chrom != s2.chrom:
        return False
    if abs(s1.pos - s2.pos) > SVConfig.CLUSTER_POS_TOL:
        return False
    if abs(s1.end - s2.end) > SVConfig.CLUSTER_POS_TOL:
        return False
    ins_len1 = int(s1.info.get('INS_LEN', 0))
    del_len1 = int(s1.info.get('DEL_LEN', 0))
    ins_len2 = int(s2.info.get('INS_LEN', 0))
    del_len2 = int(s2.info.get('DEL_LEN', 0))
    diff_ins = abs(ins_len1 - ins_len2)
    diff_del = abs(del_len1 - del_len2)
    ins_tol = SVConfig.CLUSTER_ins_LEN_TOL
    del_tol = SVConfig.CLUSTER_del_LEN_TOL
    ins_match = diff_ins <= ins_tol or diff_ins <= max(ins_len1, ins_len2) * 0.1
    del_match = diff_del <= del_tol or diff_del <= max(del_len1, del_len2) * 0.1
    if ins_match and del_match:
        return True
    return False

def _core_cluster_logic(sv_list: List[SVCandidate], match_func: Callable) -> List[SVCandidate]:
    """
    Generic linear clustering logic.
    Assumes sv_list is already sorted by (chrom, pos).
    """
    if not sv_list:
        return []
    results = []
    active = []
    for curr_sv in sorted(sv_list, key=candidate_key):
        retained = []
        for cluster in active:
            rep = cluster[0]
            if rep.chrom != curr_sv.chrom or curr_sv.pos - rep.pos > SVConfig.CLUSTER_POS_TOL:
                merged = _merge_candidates(cluster)
                if merged:
                    results.append(merged)
            else:
                retained.append(cluster)
        active = retained
        matching = [cluster for cluster in active if match_func(cluster[0], curr_sv)]
        if matching:
            min(matching, key=lambda c: candidate_key(c[0])).append(curr_sv)
        else:
            active.append([curr_sv])
    for cluster in active:
        merged = _merge_candidates(cluster)
        if merged:
            results.append(merged)
    return sorted(results, key=candidate_key)

def _merge_candidates(cluster: List[SVCandidate]) -> Optional[SVCandidate]:
    if not cluster:
        return None
    if len(cluster) == 1:
        sv = cluster[0]
        has_h1 = 'Hap1' in sv.sources
        has_h2 = 'Hap2' in sv.sources
        if has_h1 and has_h2:
            gt = '1|1'
        elif has_h1:
            gt = '1|0'
        elif has_h2:
            gt = '0|1'
        else:
            gt = './.'
        sv.genotype = gt
        sv.info['SUPPORT_H1'] = 1 if has_h1 else 0
        sv.info['SUPPORT_H2'] = 1 if has_h2 else 0
        return sv
    has_structured_id = any((sv.id for sv in cluster))
    if has_structured_id:
        candidates_with_id = [sv for sv in cluster if sv.id]
        rep = min(candidates_with_id, key=lambda sv: (len(sv.id), sv.id, candidate_key(sv)))
    else:
        cluster.sort(key=lambda x: (x.sv_len, candidate_key(x)))
        rep = cluster[len(cluster) // 2]
    sources = set()
    all_read_names = set()
    (h1_cnt, h2_cnt) = (0, 0)
    hap_names = [set(),set()]
    merged_events = set()
    merged_pvents = set()
    for sv in cluster:
        if hasattr(sv, 'read_names'):
            all_read_names.update(sv.read_names)
        if 'EVENT' in sv.info:
            merged_events.add(sv.info['EVENT'])
        if 'PEVENT' in sv.info:
            merged_pvents.add(sv.info['PEVENT'])
        sources.update(sv.sources)
        if 'Hap1' in sv.sources:
            h1_cnt += 1
            hap_names[0].update(sv.read_names)
        if 'Hap2' in sv.sources:
            h2_cnt += 1
            hap_names[1].update(sv.read_names)
    (has_h1, has_h2) = (h1_cnt > 0, h2_cnt > 0)
    if has_h1 and has_h2:
        gt = '1|1'
    elif has_h1:
        gt = '1|0'
    elif has_h2:
        gt = '0|1'
    else:
        gt = './.'
    new_sv = SVCandidate(sv_type=rep.raw_type, chrom=rep.chrom, pos=rep.pos, sv_len=rep.sv_len, read_name='', info=rep.info.copy(), end=rep.end, filter_str=rep.filter)
    new_sv.sources = sources
    new_sv.read_names = all_read_names
    new_sv.genotype = gt
    new_sv.info['SUPPORT_H1'] = len(hap_names[0]) or int(has_h1)
    new_sv.info['SUPPORT_H2'] = len(hap_names[1]) or int(has_h2)
    for key in ('ATOM_ID','CONTEXT_ID','PARENT_EVENT','REGION_ID'):
        values={part for sv in cluster for part in str(sv.info.get(key,'')).split(',') if part}
        if values:new_sv.info[key]=','.join(sorted(values))
    if any(sv.info.get('ROLE')=='SIMPLE' for sv in cluster):new_sv.info['ROLE']='SIMPLE'
    elif all(sv.info.get('ROLE')=='COMPONENT' for sv in cluster):new_sv.info['ROLE']='COMPONENT'
    new_sv.id = rep.id
    mateids = {sv.mateid for sv in cluster if sv.mateid}
    if len(mateids) == 1:
        new_sv.mateid = mateids.pop()
    elif rep.mateid:
        new_sv.mateid = rep.mateid
    if merged_events:
        new_sv.info['EVENT'] = ','.join(sorted(merged_events))
    if merged_pvents:
        new_sv.info['PEVENT'] = ','.join(sorted(merged_pvents))
    if rep.type == 'INS':
        source_fields=('TOPOLOGY','SOURCE_STATE','SOURCE_CHR','SOURCE_START','SOURCE_END','DUP_ORIENT','DUP_MODE')
        if any(sv.info.get('TOPOLOGY')=='DUP' for sv in cluster):
            proofs={tuple(sv.info.get(key) for key in source_fields) for sv in cluster}
            if len(proofs)!=1:
                # Identical inserted alleles do not imply identical retained
                # sources on both haplotypes. Keep per-haplotype graph evidence.
                for key in source_fields:new_sv.info.pop(key,None)
                new_sv.info.pop('SUBTYPE',None);new_sv.info.pop('INS_ORIGIN',None)
                new_sv.info.update(SOURCE_STATE='UNKNOWN',SOURCE_ISSUE='HAPLOTYPE_SOURCE_NOT_CONSENSUS')
        longest_seq_sv = min(cluster, key=lambda x: (-len(x.info.get('SEQ', '')), candidate_key(x)))
        if longest_seq_sv.info.get('SEQ'):
            new_sv.info['SEQ'] = longest_seq_sv.info['SEQ']
            new_sv.sv_len = longest_seq_sv.sv_len
    return new_sv

def cluster_and_genotype(hap1_svs: List[SVCandidate], hap2_svs: List[SVCandidate],report=True) -> List[SVCandidate]:
    """
    Main entry point: segregates complex events, dispatches by SV type, sorts, and clusters.
    """
    for sv in hap1_svs:
        sv.sources = {'Hap1'}
        sv.genotype = '1|0'
    for sv in hap2_svs:
        sv.sources = {'Hap2'}
        sv.genotype = '0|1'
    raw_svs = hap1_svs + hap2_svs
    if not raw_svs:
        return []
    final_results = []
    clusterable_buckets = defaultdict(list)
    for sv in raw_svs:
        if sv.filter == 'COMPLEX':
            final_results.append(sv)
        elif sv.type=='CPX':
            clusterable_buckets['CPX'].append(sv)
        else:
            clusterable_buckets[sv.type].append(sv)
    common_sort_key = candidate_key
    for sv_type in ['DEL', 'DUP', 'INV']:
        if (candidates := clusterable_buckets[sv_type]):
            candidates.sort(key=common_sort_key)
            final_results.extend(_core_cluster_logic(candidates, _match_span_sv))
    if (ins_candidates := clusterable_buckets['INS']):
        ins_candidates.sort(key=common_sort_key)
        final_results.extend(_core_cluster_logic(ins_candidates, _match_ins_sv))
    for sv_type in ['BND']:
        if (bnd_tra_candidates := clusterable_buckets[sv_type]):
            bnd_tra_candidates.sort(key=lambda x: (x.chrom, x.pos, x.info.get('CHR2', ''), x.info.get('POS2', 0)))
            final_results.extend(_core_cluster_logic(bnd_tra_candidates, _match_bnd_tra_sv))
    if clusterable_buckets.get('DELINS'):
        final_results.extend(_core_cluster_logic(clusterable_buckets['DELINS'],_match_composite_sv))
    if clusterable_buckets.get('CPX'):
        final_results.extend(_core_cluster_logic(clusterable_buckets['CPX'],_match_composite_sv))
    if report:logger.info(f'Genotyping Summary: Input={len(raw_svs)}, Output={len(final_results)}')
    final_results.sort(key=candidate_key)
    return final_results


def cluster_structural_calls(hap1_svs,hap2_svs):
    """Keep copy-local components out of reference-locus event genotypes.

    The original INS/DEL match functions are unchanged. Components are compared
    only inside the same confirmed parent and coordinate/occurrence context.
    Reference-indel projections are independently clustered by the caller.
    """
    primary=cluster_and_genotype([c for c in hap1_svs if c.info.get('ROLE') not in ('COMPONENT','PROJECTION')],
        [c for c in hap2_svs if c.info.get('ROLE') not in ('COMPONENT','PROJECTION')])
    aliases={alias:c.id for c in primary if c.id for alias in str(c.info.get('EVENT','')).split(',') if alias}
    groups=defaultdict(lambda:[[],[]])
    for hap,calls in enumerate((hap1_svs,hap2_svs)):
        for c in calls:
            if c.info.get('ROLE')!='COMPONENT':continue
            parents=','.join(sorted({aliases.get(p,p) for p in c.info['PARENT_EVENT'].split(',')}))
            c.info['PARENT_EVENT']=parents
            groups[(parents,c.info.get('ALLELE_CONTEXT','REFERENCE'))][hap].append(c)
    result=list(primary)
    for key,haps in sorted(groups.items()):result.extend(cluster_and_genotype(*haps,report=False))
    logger.info('Scoped component genotyping: groups=%d output=%d',len(groups),len(result)-len(primary))
    for c in primary:
        if c.info.get('ROLE')=='EVENT':c.info['EVENT']=c.id
    return sorted(result,key=candidate_key),aliases
