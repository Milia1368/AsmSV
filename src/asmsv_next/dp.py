"""An additive global path objective followed by backbone selection inside that path."""
from dataclasses import dataclass
import math
from asmsv_next.config import SVConfig
from asmsv_next.constants import (DP_TERMINAL_RATE, DP_EXCESS_OVERLAP_START,
    DP_REFERENCE_RATE, DP_REFERENCE_CAP, DP_WEAK_TRA_PREMIUM_RATIO)


def node_key(b):
    return (b.qs,b.qe,b.ref_name,b.strand,b.ref_start,b.ref_end,tuple(b.cigar),b.as_score)


def oriented(b):
    return (b.ref_start,b.ref_end) if b.strand=='+' else (-b.ref_end,-b.ref_start)


def legal(a,b,collinear=False):
    if not (a.qs<b.qs and a.qe<b.qe): return False
    if not collinear:return True
    u,v=oriented(a);x,y=oriented(b)
    return a.ref_name==b.ref_name and a.strand==b.strand and u<x and v<y


def chromosome_cost(a,b,p):
    base=p.cross_chromosome
    if not base:return 0.
    premium=base*DP_WEAK_TRA_PREMIUM_RATIO
    # One and two boundary budgets, measured in AQ information rather than bp.
    return base+premium*min(1.,max(0.,2-min(a.as_score,b.as_score)/base))


def edge_terms(a,b,preset=None):
    p=preset or SVConfig.CALLER_PRESET.dp
    gq=b.qs-a.qe;overlap=max(0,-gq)
    shared=max(0.0,overlap_reward(a,b.qs,a.qe),overlap_reward(b,b.qs,a.qe)) if overlap else 0.0
    terms=dict(query_gap=p.query_gap*max(0,gq),shared_query=shared,
        excess_overlap=p.excess_overlap_rate*max(0,overlap-DP_EXCESS_OVERLAP_START),
        chromosome_change=0.0,orientation_change=0.0,reference_residual=0.0)
    # Charge every observed boundary; chromosome and orientation can coexist.
    if a.ref_name!=b.ref_name:terms['chromosome_change']=chromosome_cost(a,b,p)
    if a.strand!=b.strand:
        terms['orientation_change']=p.inversion
        if a.ref_name==b.ref_name:
            exit_a=a.ref_end if a.strand=='+' else a.ref_start
            entry_b=b.ref_start if b.strand=='+' else b.ref_end
            displacement=max(0,abs(entry_b-exit_a)-max(a.ref_end-a.ref_start,b.ref_end-b.ref_start))
            terms['reference_residual']=min(DP_REFERENCE_CAP,DP_REFERENCE_RATE*displacement)
    elif a.ref_name==b.ref_name:
        gr=oriented(b)[0]-oriented(a)[1]
        terms['reference_residual']=min(DP_REFERENCE_CAP,DP_REFERENCE_RATE*abs(gq-gr))
    return terms


def overlap_reward(block,start,end):
    """Signed local AQ reward; floor rescues have only their small uniform credit.

    NM-only evidence remains an explicitly estimated slice, never a positional
    match claim. Simple external nodes fall back to their declared density.
    """
    root=getattr(block,'score_root',None)
    if root is None or getattr(block,'reward_rescued',False):
        return (end-start)*block.as_score/(block.qe-block.qs)
    return root.index.score(start,end).node_reward*block.mapq_factor


def transition_cost(a,b,preset=None): return math.fsum(edge_terms(a,b,preset).values())


def path_score(blocks,query_length,preset=None):
    p=preset or SVConfig.CALLER_PRESET.dp
    if not blocks:return -DP_TERMINAL_RATE*query_length
    return math.fsum([b.as_score for b in blocks]+
        [-transition_cost(a,b,p) for a,b in zip(blocks,blocks[1:])]+
        [-DP_TERMINAL_RATE*blocks[0].qs,-DP_TERMINAL_RATE*(query_length-blocks[-1].qe)])


@dataclass(frozen=True)
class PathResult:
    blocks: tuple
    score: float
    evaluated_predecessors: int
    tied: bool


def solve_path(nodes,query_length,*,collinear=False,preset=None):
    p=preset or SVConfig.CALLER_PRESET.dp
    if type(query_length) is not int or query_length<0:raise ValueError('Invalid contig length')
    blocks=sorted((b for b in nodes if b.is_valid and b.as_score>0),key=node_key)
    if any(not 0<=b.qs<b.qe<=query_length or not math.isfinite(b.as_score) for b in blocks):
        raise ValueError('Invalid DP node coordinates/reward')
    unique={}
    for b in blocks:unique.setdefault(node_key(b),b)
    blocks=list(unique.values());scores=[];parents=[];sizes=[];counts=[];evaluated=0
    for j,b in enumerate(blocks):
        best=b.as_score-DP_TERMINAL_RATE*b.qs;parent=-1;size=1;count=1
        for i in range(j):
            a=blocks[i]
            if not legal(a,b,collinear):continue
            evaluated+=1
            value=scores[i]+b.as_score-transition_cost(a,b,p)
            if value>best:count=counts[i]
            elif value==best:count=min(2,count+counts[i])
            if value>best or (value==best and (sizes[i]+1,i)<(size,parent)):
                best,parent,size=value,i,sizes[i]+1
        scores.append(best);parents.append(parent);sizes.append(size);counts.append(count)
    best=-DP_TERMINAL_RATE*query_length;end=-1;count=1
    for i,b in enumerate(blocks):
        value=scores[i]-DP_TERMINAL_RATE*(query_length-b.qe)
        if value>best:count=counts[i]
        elif value==best:count=min(2,count+counts[i])
        if value>best or (value==best and end>=0 and (sizes[i],i)<(sizes[end],end)):
            best,end=value,i
    chosen=[]
    while end>=0:chosen.append(blocks[end]);end=parents[end]
    chosen=tuple(reversed(chosen));actual=path_score(chosen,query_length,p)
    if not math.isclose(actual,best,rel_tol=1e-10,abs_tol=1e-7):
        raise ArithmeticError('DP traceback and objective disagree')
    return PathResult(chosen,actual,evaluated,count>1)


def select_backbone(chain,query_length):
    for b in chain:b.is_backbone=False
    result=solve_path(chain,query_length,collinear=True)
    selected={id(b) for b in result.blocks}
    indices={i for i,b in enumerate(chain) if id(b) in selected}
    for i in indices:chain[i].is_backbone=True
    return indices,result.blocks[0].strand if result.blocks else '+'
