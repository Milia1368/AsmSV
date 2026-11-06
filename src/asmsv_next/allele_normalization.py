"""Exact allele-equivalence transformations; no truth coordinates are used."""
from copy import copy
import hashlib
from asmsv_next.sequence_check import reverse_complement


def clone_candidate(candidate):
    result=copy(candidate);result.info=candidate.info.copy()
    result.sources=set(candidate.sources);result.read_names=set(candidate.read_names)
    result.format_data=candidate.format_data.copy()
    return result


def shifted_insertion(reference,chrom,position,sequence,target):
    """Return the rotated payload only when both complete local alleles agree."""
    if target==position:return sequence
    if min(position,target)<0:return None
    try:pad=reference.fetch(chrom,min(position,target),max(position,target)).upper()
    except (ValueError,KeyError):return None
    if len(pad)!=abs(target-position) or any(c not in 'ACGT' for c in pad+sequence):return None
    if target>position:
        allele=sequence+pad
        return allele[len(pad):] if allele.startswith(pad) else None
    allele=pad+sequence
    return allele[:-len(pad)] if allele.endswith(pad) else None


def canonical_tandem_indels(calls,parents,reference):
    """Normalize only the reference-change projection; original atoms stay intact."""
    answer=[]
    for candidate in calls:
        if candidate.type!='INS' or not candidate.info.get('SEQ'):
            answer.append(candidate);continue
        matches=[parents[p] for p in str(candidate.info.get('PARENT_EVENT','')).split(',')
                 if p in parents and parents[p].type=='DUP' and parents[p].info.get('DUP_MODE')=='TANDEM']
        if len(matches)!=1:
            answer.append(candidate);continue
        parent=matches[0];target=parent.info.get('TARGET_POS')
        if parent.info.get('TARGET_CHR')!=candidate.chrom or target is None or target==candidate.pos:
            answer.append(candidate);continue
        original=candidate.info['SEQ'];rotated=shifted_insertion(reference,candidate.chrom,candidate.pos,original,int(target))
        if rotated is None:
            answer.append(candidate);continue
        result=clone_candidate(candidate)
        result.pos=result.end=int(target);result.info['SEQ']=rotated
        result.info.update(ALLELE_NORMALIZATION='EXACT_TANDEM_TARGET',
            ORIGINAL_POS=candidate.pos,ORIGINAL_SEQUENCE_SHA256=hashlib.sha256(original.encode()).hexdigest())
        answer.append(result)
    return answer


def canonical_insertions(calls,reference,*,indel_view=False,limit=15):
    """Bounded left normalization with a certified natural stop.

    Applies to pure ACGT insertions only. Long repeat boundaries extending
    beyond the budget are left unchanged. Coordinates of raw graph atoms and
    physical query spans are never rewritten by this output transformation.
    """
    answer=[]
    for call in calls:
        sequence=call.info.get('SEQ','')
        if (call.type!='INS' or not sequence or call.end!=call.pos or
                not indel_view and call.info.get('ROLE') in ('COMPONENT','PROJECTION') or
                any(c not in 'ACGT' for c in sequence)):
            answer.append(call);continue
        p=call.pos
        try:flank=reference.fetch(call.chrom,max(1,p-limit-1),p).upper()
        except (ValueError,KeyError):answer.append(call);continue
        delta=0;n=len(sequence)
        while delta<len(flank) and flank[-delta-1]==sequence[(n-delta-1)%n]:delta+=1
        if not delta or delta>limit:answer.append(call);continue
        rotated=shifted_insertion(reference,call.chrom,p,sequence,p-delta)
        if rotated is None:answer.append(call);continue
        result=clone_candidate(call);result.pos=result.end=p-delta
        result.info.update(SEQ=rotated,ALLELE_NORMALIZATION='EXACT_BOUNDED_LEFT',
            ORIGINAL_POS=p,ORIGINAL_SEQUENCE_SHA256=hashlib.sha256(sequence.encode()).hexdigest())
        answer.append(result)
    return answer


def canonical_inversion(reference,chrom,start,end,sequence,qspan,strand,limit=15):
    """Extend an inversion through a bounded exact complementary boundary.

    Require a natural stop within the limit and actual query support for both
    added flanks. This yields the same maximal microhomologous interval for
    split-block and sequence-derived evidence without fuzzy genotype merging.
    """
    if reference is None or end<=start:return None
    try:
        length=reference.get_reference_length(chrom)
        left=reference.fetch(chrom,max(0,start-limit-1),start).upper()
        right=reference.fetch(chrom,end,min(length,end+limit+1)).upper()
    except (ValueError,KeyError):return None
    delta=0
    while delta<min(len(left),len(right)) and left[-delta-1] in 'ACGT' and right[delta] in 'ACGT' and left[-delta-1]==reverse_complement(right[delta]):
        delta+=1
    if not delta or delta>limit:return None
    qs,qe=qspan
    if qs<delta or qe+delta>len(sequence):return None
    a=left[-delta:];b=right[:delta]
    prefix,suffix=(a,b) if strand=='+' else (reverse_complement(b),reverse_complement(a))
    if sequence[qs-delta:qs].upper()!=prefix or sequence[qe:qe+delta].upper()!=suffix:return None
    return start-delta,end+delta,(qs-delta,qe+delta)
