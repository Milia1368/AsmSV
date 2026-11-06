"""Bounded exact sequence recovery between already selected anchors.

No aligner is imported on the default path. Coordinates returned here are local,
0-based interbase coordinates along the query-oriented reference residual.
"""
from collections import defaultdict
from dataclasses import dataclass, field
from fractions import Fraction
from functools import reduce
from math import gcd, lcm

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
MAX_LINEAR_BASES = 2000000
SPARSE_MAX_POINTS = 20000
SPARSE_MASK = 127


def reverse_complement(sequence):
    return sequence.translate(RC)[::-1]


@dataclass
class SequenceCheck:
    status: str
    changes: list = field(default_factory=list)
    homologies: list = field(default_factory=list)
    issues: list = field(default_factory=list)
    method: str = 'EXACT'
    alignment: list = field(default_factory=list)
    alternative_homologies: list = field(default_factory=list)


def shared_seed(query, reference, k):
    if len(query) > len(reference):query,reference=reference,query
    seeds={query[i:i+k] for i in range(len(query)-k+1)}
    return any(reference[i:i+k] in seeds for i in range(len(reference)-k+1))


def unique_chain(query, reference, k=16, max_blocks=1000):
    """Unique exact seeds -> diagonal runs -> maximum covered collinear chain.

    A transition trims the next exact run equally on both axes. This permits
    boundary microhomology without assigning a query base to two homologies.
    """
    def positions(sequence):
        result={}
        for i in range(len(sequence)-k+1):
            seed=sequence[i:i+k]
            result[seed]=None if seed in result else i
        return result
    qp,rp=positions(query),positions(reference)
    points=sorted((q,rp[seed]) for seed,q in qp.items()
                  if q is not None and rp.get(seed) is not None)
    blocks=[]
    for q,r in points:
        if blocks and q-r==blocks[-1][0]-blocks[-1][2] and q<=blocks[-1][1]:
            blocks[-1][1]=max(blocks[-1][1],q+k)
            blocks[-1][3]=max(blocks[-1][3],r+k)
        else:blocks.append([q,q+k,r,r+k])
    return _chain_blocks(blocks,k,max_blocks)


def _chain_blocks(blocks,k,max_blocks):
    if len(blocks)>max_blocks:return [],'SEED_BUDGET'
    scores=[];parents=[];trims=[]
    for j,b in enumerate(blocks):
        best=b[1]-b[0];parent=-1;trim=0
        for i,a in enumerate(blocks[:j]):
            overlap=max(0,a[1]-b[0],a[3]-b[2]);remaining=b[1]-b[0]-overlap
            if b[1]>a[1] and b[3]>a[3] and remaining>=k and scores[i]+remaining>best:
                best=scores[i]+remaining;parent=i;trim=overlap
        scores.append(best);parents.append(parent);trims.append(trim)
    index=max(range(len(scores)),key=scores.__getitem__) if scores else -1
    chain=[];used=set()
    while index>=0:
        b=list(blocks[index]);b[0]+=trims[index];b[2]+=trims[index]
        chain.append(b);used.add(index);index=parents[index]
    chain.reverse()
    return chain,('SEED_ORDER_CONFLICT' if len(used)<len(blocks) else 'OK')


def sparse_unique_chain(query,reference,k=16,max_blocks=1000,*,return_blocks=False):
    """Content-selected exact seeds with bounded scan, points, and block DP.

    Selection depends on seed content, so all copies of a selected seed are
    counted even across shifted strings. A selected seed must occur once in
    each string. Maximal exact extensions, rather than sampled bases alone,
    establish the restored homology. Limits cause abstention.
    """
    if k!=16 or len(query)+len(reference)>MAX_LINEAR_BASES:return [],'SEED_BUDGET'
    def positions(sequence):
        result={};code=0;mask=(1<<(2*k))-1;alphabet={'A':0,'C':1,'G':2,'T':3};count=0
        for i,base in enumerate(sequence):
            code=((code<<2)|alphabet[base])&mask
            if i<k-1:continue
            mixed=code
            mixed=((mixed^(mixed>>16))*0x7feb352d)&0xffffffff
            mixed=((mixed^(mixed>>15))*0x846ca68b)&0xffffffff
            mixed^=mixed>>16
            if mixed&SPARSE_MASK:continue
            count+=1
            if count>SPARSE_MAX_POINTS:return None
            result[code]=None if code in result else i-k+1
        return result
    qp,rp=positions(query),positions(reference)
    if qp is None or rp is None:return [],'SEED_BUDGET'
    diagonals=defaultdict(list)
    for seed,q in qp.items():
        r=rp.get(seed)
        if q is not None and r is not None:diagonals[q-r].append((q,r))
    blocks=[]
    for points in diagonals.values():
        covered=-1
        for q,r in sorted(points):
            if q<covered:continue
            a,b=q,r;s,t=q+k,r+k
            while a and b and query[a-1]==reference[b-1]:a-=1;b-=1
            while s<len(query) and t<len(reference) and query[s]==reference[t]:s+=1;t+=1
            blocks.append([a,s,b,t]);covered=s
            if len(blocks)>max_blocks:return [],'SEED_BUDGET'
    if not blocks:return [],'SEED_BUDGET'
    chain,status=_chain_blocks(sorted(blocks),k,max_blocks)
    return (chain,status,sorted(blocks)) if return_blocks else (chain,status)


def exact_filler(query, reference, side, k=16):
    """A complete filler can establish target restoration at one reference edge."""
    if len(query)<k or any(c not in 'ACGT' for c in query+reference):return False
    if not (reference.startswith(query) if side=='right' else reference.endswith(query)):return False
    first=reference.find(query)
    return first>=0 and reference.find(query,first+1)<0


def common_prefix(a,b):
    i=0;limit=min(len(a),len(b))
    while i<limit and a[i]==b[i]:i+=1
    return i


def exact_flank_model(query,reference,*,k=16):
    """Exact homology, pure indel or inversion with unchanged forward flanks.

    The complete compared allele is explained. At most k-1 bases of boundary
    microhomology are reconsidered; no seed DP or similarity threshold is used.
    Returned coordinates refer to the original strings, including any padding.
    """
    if query.strip('ACGT') or reference.strip('ACGT'):return None
    if query==reference:
        return SequenceCheck('HOMOLOGY_RESTORED',homologies=[[0,len(query),0,len(reference)]] if query else [])
    p=common_prefix(query,reference)
    s=common_prefix(query[p:][::-1],reference[p:][::-1])
    def homologies(left,right):
        blocks=[]
        if left:blocks.append([0,left,0,left])
        if right:blocks.append([len(query)-right,len(query),len(reference)-right,len(reference)])
        return blocks
    qe,re=len(query)-s,len(reference)-s
    if p==qe or p==re:
        return SequenceCheck('EXPLAINED',changes=[dict(type='GAP',query=[p,qe],reference=[p,re])],
            homologies=homologies(p,s),method='EXACT_FLANKS')
    if len(query)!=len(reference):return None
    reverse=reverse_complement(reference);models=[]
    # Every admissible inversion contains the same untrimmed central core.
    # Locate that core once in the bounded reverse-reference window, then
    # compare only the at-most k-1 flanking bases for each boundary choice.
    # This avoids copying/comparing an entire long core for every k*k pair.
    core=query[p:qe];minimum_right=max(0,s-k+1)
    hit=reverse.find(core,minimum_right,min(len(reverse),s+k-1+len(core)))
    while hit>=0:
        for left in range(max(0,p-k+1),p+1):
            right=hit-p+left;end=len(query)-right
            if (minimum_right<=right<=s and end-left>=k and
                    query[left:p]==reverse[right:hit] and
                    query[qe:end]==reverse[hit+len(core):len(reference)-left]):
                models.append((left,right))
        hit=reverse.find(core,hit+1,min(len(reverse),s+k-1+len(core)))
    if not models:return None
    left,right=min(models,key=lambda pair:(sum(pair),pair))
    end=len(query)-right
    return SequenceCheck('EXACT_INV' if not left and not right else 'EXPLAINED',
        changes=[dict(type='INV',query=[left,end],reference=[left,end])],
        homologies=homologies(left,right),method='EXACT_FLANKS')


def exact_inverse_edit(query,reference,*,k=16,min_edit=30,max_bases=40000,min_coverage=.8):
    """One internal SV edit in an otherwise exact inverted interval.

    This bounded linear template does not run a seed DP or an aligner. Forward
    identical flanks are removed, with <k bp boundary wobble reconsidered.
    Both inverse arms must be long and unique within the compared strings.
    """
    if len(query)+len(reference)>max_bases:return None
    query,reference=query.upper(),reference.upper()
    if any(c not in 'ACGT' for c in query+reference):return None
    p=common_prefix(query,reference)
    s=common_prefix(query[p:][::-1],reference[p:][::-1])
    reverse=reverse_complement(reference);arm_min=max(2*k,min_edit);models=[]
    for left in range(max(0,p-k+1),p+1):
        for right in range(max(0,s-k+1),s+1):
            q=query[left:len(query)-right if right else len(query)]
            t=reverse[right:len(reference)-left]
            if min(len(q),len(t))<2*arm_min or q[:arm_min]!=t[:arm_min] or q[-arm_min:]!=t[-arm_min:]:continue
            if q==t:continue
            a=common_prefix(q,t);b=common_prefix(q[a:][::-1],t[a:][::-1])
            if min(a,b)<arm_min or (a+b)/max(len(q),len(t))<min_coverage:continue
            qe,te=len(q)-b,len(t)-b;ins,deleted=qe-a,te-a
            if not ((ins==0 and deleted>=min_edit) or (deleted==0 and ins>=min_edit) or min(ins,deleted)>=min_edit):continue
            if ins and deleted and shared_seed(q[a:qe],t[a:te],k):continue
            arms=(q[:a],q[qe:])
            if any(text.find(arm,text.find(arm)+1)>=0 for arm in arms for text in (q,t)):continue
            f=reference[left:len(reference)-right if right else len(reference)]
            x=common_prefix(q,f);y=common_prefix(q[x:][::-1],f[x:][::-1])
            if min(x,y)>=arm_min and (x+y)/max(len(q),len(f))>=min_coverage:continue
            models.append(dict(prefix=left,suffix=right,query=[left,len(query)-right],
                reference=[left,len(reference)-right],query_edit=[left+a,left+qe],
                reference_edit=[left+len(t)-te,left+len(t)-a],
                inverse_arms=[[left,left+a,left+len(t)-a,left+len(t)],
                              [left+qe,len(query)-right,left,left+b]],
                ins_len=ins,del_len=deleted,homologous_fraction=(a+b)/max(len(q),len(t))))
    if not models:return None
    # Maximal equivalent interval is a deterministic representation, not a
    # claim that microhomologous historical breakpoints are uniquely known.
    return min(models,key=lambda m:(-(m['reference'][1]-m['reference'][0]),m['prefix'],m['suffix'],m['reference_edit']))


def check_sequence(query, reference, *, k=16, max_bases=20000, max_blocks=1000,
                   min_gap=100, max_ref_gap=150000, realign='off', preset=None):
    query,reference=query.upper(),reference.upper()
    if any(c not in 'ACGT' for c in query):return SequenceCheck('ASSEMBLY_GAP',issues=['ASSEMBLY_GAP'])
    if any(c not in 'ACGT' for c in reference):return SequenceCheck('REFERENCE_GAP',issues=['REFERENCE_GAP'])
    if query==reference:
        return SequenceCheck('HOMOLOGY_RESTORED',homologies=[[0,len(query),0,len(reference)]])
    if len(query)>=k and query==reverse_complement(reference):
        return SequenceCheck('EXACT_INV',changes=[dict(type='INV',query=[0,len(query)],reference=[0,len(reference)])])
    sparse=len(query)+len(reference)>max_bases
    if sparse:
        exact=exact_flank_model(query,reference,k=k)
        if exact is not None:return exact
    alternatives=[]
    if sparse:
        answer=sparse_unique_chain(query,reference,k,max_blocks,return_blocks=True)
        chain,status=answer[:2]
        if len(answer)==3 and status=='SEED_ORDER_CONFLICT':
            # Short off-path repeat matches cannot veto a much larger exact
            # island when every remaining edit is an exact pure indel.
            alternatives=[b for b in answer[2] if not any(c[0]<=b[0] and b[1]<=c[1] and c[2]<=b[2] and b[3]<=c[3] and c[0]-c[2]==b[0]-b[2] for c in chain)]
            pure=bool(chain);qs=rs=0
            for block in chain+[[len(query),len(query),len(reference),len(reference)]]:
                q,r=query[qs:block[0]],reference[rs:block[2]]
                p=common_prefix(q,r);s=common_prefix(q[p:][::-1],r[p:][::-1])
                if len(q)-p-s and len(r)-p-s:pure=False
                qs,rs=block[1],block[3]
            if (pure and alternatives and all(b[1]-b[0]<min_gap for b in alternatives) and
                    sum(b[1]-b[0] for b in alternatives)<max(b[1]-b[0] for b in chain)):
                status='OK'
    else:chain,status=unique_chain(query,reference,k,max_blocks)
    if status!='OK':
        result=SequenceCheck(status,homologies=chain,issues=[status],method='KMER')
    else:
        changes=[];issues=[];qs=rs=0
        for block in chain+[[len(query),len(query),len(reference),len(reference)]]:
            qe,re=block[0],block[2]
            if qe>qs or re>rs:
                q,r=query[qs:qe],reference[rs:re]
                if q==r:pass
                elif len(q)>=k and q==reverse_complement(r):
                    changes.append(dict(type='INV',query=[qs,qe],reference=[rs,re]))
                elif len(q)>=min_gap and len(r)>=min_gap:
                    if sparse and len(q)+len(r)>max_bases:issues.append('LARGE_RESIDUAL_BUDGET')
                    elif len(r)<=max_ref_gap and not shared_seed(q,r,k):
                        changes.append(dict(type='DELINS',query=[qs,qe],reference=[rs,re],sequence=q))
                    else:issues.append('SEED_SHARED_BUT_UNRESOLVED')
                else:
                    changes.append(dict(type='GAP',query=[qs,qe],reference=[rs,re]))
            qs,rs=block[1],block[3]
        result=SequenceCheck('PARTIAL' if issues else 'EXPLAINED',changes,chain,issues,'KMER')
    if sparse:result.method='SPARSE_EXACT_KMER'
    result.alternative_homologies=alternatives
    if result.issues and result.status!='SEED_ORDER_CONFLICT' and realign=='fallback':
        return _realign(query,reference,preset,k,max_bases,result)
    return result


def _realign(query, reference, preset, k, max_bases, previous):
    """Optional derived evidence; never reprocesses an original alignment view."""
    if len(query)+len(reference)>max_bases:return previous
    try:
        from pywfa import WavefrontAligner
    except ImportError:
        previous.issues.append('REALIGN_UNAVAILABLE');return previous
    # Native match reward 1 is removed using the fixed-length identity
    # 2M=|Q|+|R|-2X-sum(gap lengths). No SSV acceptance threshold changes.
    p=preset.aq
    costs=[Fraction(str(2*(1+p.mismatch)))]+[Fraction(str(2*x)) for x in p.gap_open]+[Fraction(str(2*x+1)) for x in p.gap_extend]
    scale=reduce(lcm,(v.denominator for v in costs),1)
    integers=[int(v*scale) for v in costs];common=reduce(gcd,integers)
    integers=[v//common for v in integers]
    mismatch=integers[0];n=len(p.gap_open)
    opens=integers[1:1+n];extends=integers[1+n:]
    if n!=2:
        previous.issues.append('REALIGN_MODEL_UNSUPPORTED');return previous
    aligner=WavefrontAligner(reference,span='end-to-end',distance='affine2p',match=0,
        mismatch=mismatch,gap_opening=opens[0],gap_extension=extends[0],
        gap_opening2=opens[1],gap_extension2=extends[1],heuristic=None,
        memory_mode='biwfa',max_steps=200000,max_memory=128*1024*1024)
    answer=aligner(query,clip_cigar=False)
    if answer.status!=0:
        previous.issues.append('RECHECK_LIMIT');return previous
    from asmsv_next.alignment import AlignedSegment
    cigar=[(7,64)]+[(7 if op==0 else op,length) for op,length in answer.cigartuples]+[(7,64)]
    nm=sum(length for op,length in cigar if op in (1,2,8))
    block=AlignedSegment(dict(q_name='recheck',ref_name='local',ref_start=0,
        ref_end=len(reference)+128,cigar=cigar,mapq=60,query_length=len(query)+128,
        tags={'NM':nm},source={'aligner':preset.aligner,'mode':preset.mode}))
    from asmsv_next.config import SVConfig
    atoms=block.parse_intra_alignment_sv('A'*64+query+'C'*64,
        SVConfig.MIN_SV_SIZE,SVConfig.MERGE_INDEL_TOL,with_provenance=True)
    changes=[]
    for atom in atoms:
        parts=atom['evidence']['parts'];q0=min(x[1] for x in parts)-64;q1=max(x[2] for x in parts)-64
        r0=atom['pos']-64;r1=r0+(atom['len'] if atom['type']=='DEL' else 0)
        changes.append(dict(type=atom['type'],query=[q0,q1],reference=[r0,r1],
            length=atom['len'],sequence=atom['info'].get('SEQ',''),derived_cigar=True))
    return SequenceCheck('EXPLAINED',changes,[],[],'WFA',answer.cigartuples)
