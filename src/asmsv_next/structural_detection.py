"""Finite, structure-first interpretation of a selected chain.

Original SSV kernels generate atoms. This module owns their scopes, scans each
selected CIGAR view once, and assigns a simple event, a supported CSV, or an
explicitly unresolved context. It does not infer mutation chronology.
"""
from bisect import bisect_left,bisect_right
from collections import defaultdict
from dataclasses import dataclass, field
import hashlib
import json
import time
from types import SimpleNamespace

from asmsv_next.config import SVConfig
from asmsv_next.sequence_check import check_sequence,exact_filler,reverse_complement,shared_seed,exact_inverse_edit,common_prefix,exact_flank_model
from asmsv_next.ssv import make_gap_caller
from asmsv_next.structures import SVCandidate
from asmsv_next.nomenclature import duplication_metadata,complex_metadata
from asmsv_next.sequence_check import MAX_LINEAR_BASES


def identifier(prefix,*parts):
    payload=json.dumps(parts,sort_keys=True,separators=(',',':'),ensure_ascii=False)
    return prefix+'.'+hashlib.sha256(payload.encode()).hexdigest()[:24]


def junction(haplotype,name,context,left,right,query=None,atom_id=None):
    ports=sorted((left,right),key=lambda p:(p['chrom'],p['position'],p['side']))
    return dict(junction_id=identifier('JUNC',haplotype,name,context,ports,atom_id),left=left,right=right,
        query_support=query,atom_id=atom_id)


def union_length(intervals):
    end=None;total=0
    for a,b in sorted(intervals):
        if b<=a:continue
        if end is None or a>end:total+=b-a
        elif b>end:total+=b-end
        end=max(end if end is not None else b,b)
    return total


def mapped_parts(block,query=None,reference=None,exact=False):
    """Read the existing AQ index, never rerun the CIGAR SV state machine."""
    q0,q1=query or (block.qs,block.qe)
    q0,q1=max(q0,block.qs),min(q1,block.qe)
    if q1<=q0:return
    for run in block.score_root.index.local_runs(q0,q1):
        if run.op not in ((7,) if exact else (0,7,8)):continue
        s,e=max(q0,run.qs),min(q1,run.qe)
        if s>=e:continue
        if block.strand=='+':a,b=run.rs+s-run.qs,run.rs+e-run.qs
        else:a,b=run.re-(e-run.qs),run.re-(s-run.qs)
        if reference is not None:
            lo,hi=max(a,reference[0]),min(b,reference[1])
            if lo>=hi:continue
            if block.strand=='+':s,e=s+lo-a,e-(b-hi)
            else:s,e=s+b-hi,e-(lo-a)
            a,b=lo,hi
        yield s,e,a,b


def map_port(block,query_point):
    if not block.qs<=query_point<=block.qe:return None
    index=block.score_root.index
    # Include the left aligned run at an exact boundary, even if one or more
    # zero-query deletions share that coordinate. This preserves the original
    # first-in-query-order endpoint convention.
    lo=bisect_left(index.ends,query_point);hi=bisect_right(index.starts,query_point)
    for run in index.runs[lo:hi]:
        if run.op not in (0,7,8):continue
        s,e=max(block.qs,run.qs),min(block.qe,run.qe)
        if s>=e:continue
        if s<=query_point<=e:
            return run.rs+query_point-run.qs if block.strand=='+' else run.re-(query_point-run.qs)
    return None


def entry(block):return block.ref_start if block.strand=='+' else block.ref_end
def exit_port(block):return block.ref_end if block.strand=='+' else block.ref_start


@dataclass(slots=True)
class Atom:
    atom_id: str
    candidate: object
    query: tuple
    reference_parts: list
    parts: list
    block_index: object = None
    contexts: set = field(default_factory=set)
    parent: str = ''
    explained: bool = True
    source_query: object = None


@dataclass(slots=True)
class DetectionResult:
    calls: list
    contexts: list
    events: list
    atoms: list
    cigar_scans: int
    cigar_seconds: float
    projections: list = field(default_factory=list)
    regions: list = field(default_factory=list)
    region_seconds: float = 0.


class Detector:
    def __init__(self,chain,strand,backbone,read_name,sequence,reference=None,haplotype=''):
        self.chain=chain;self.strand=strand;self.backbone=sorted(backbone)
        self.name=read_name;self.sequence=sequence;self.reference=reference;self.haplotype=haplotype
        self.kernel=make_gap_caller(read_name,sequence)
        self.cache={};self.atoms={};self.gap_cache={}
        self.parents=[];self.contexts=[];self.events=[];self.cigar_seconds=0.
        self.projections=[]
        self.block_contexts=defaultdict(list)
        self.k=SVConfig.CALLER_PRESET.aq.min_anchor
        self.source_index=defaultdict(list)
        for i in self.backbone:
            b=chain[i]
            for q0,q1,r0,r1 in mapped_parts(b):self.source_index[b.ref_name].append((r0,r1,q0,q1,i))
        self.source_starts={};self.source_max_ends={}
        for chrom,parts in self.source_index.items():
            parts.sort();self.source_starts[chrom]=[p[0] for p in parts]
            maximum=0;ends=[]
            for part in parts:maximum=max(maximum,part[1]);ends.append(maximum)
            self.source_max_ends[chrom]=ends

    def refseq(self,chrom,start,end,strand='+'):
        if self.reference is None or end<start:return None
        try:sequence=self.reference.fetch(chrom,max(0,start),end).upper()
        except (ValueError,KeyError):return None
        if len(sequence)!=end-start:return None
        return sequence if strand=='+' else reverse_complement(sequence)

    def context(self,left,right,members,kind='UNKNOWN',scope=None):
        parts=(self.haplotype,self.name,left,right,members)
        key=identifier('CTX',*parts) if scope is None else identifier('CTX',*parts,scope)
        context=dict(context_id=key,query_name=self.name,haplotype=self.haplotype,
            left_backbone=left,right_backbone=right,member_block_indices=members,
            frame=kind,issues=[],separators=[],atom_ids=[],resolution='PARTIAL',
            decomposition='UNDETERMINED',evidence_scope='SELECTED_PATH',
            alt_assessment='NOT_AVAILABLE',region_annotations=[])
        if left is None:context['issues'].append('OPEN_LEFT')
        if right is None:context['issues'].append('OPEN_RIGHT')
        if scope is not None:context['scope_atom_id']=scope
        self.contexts.append(context)
        for index in members:
            self.block_contexts[index].append(context)
            for key in self.cache.get(index,[]):self.annotate_context(self.atoms[key],context)
        return context

    def atom(self,candidate,query,reference_parts,parts,block_index=None,context=None,explained=True):
        root=self.chain[block_index].score_root.raw.alignment_id if block_index is not None else None
        key=identifier('ATOM',self.haplotype,self.name,root,candidate.type,parts)
        if key not in self.atoms:
            candidate.info.update(ATOM_ID=key,ROLE='SIMPLE',EVIDENCE_SCOPE='SELECTED_PATH')
            self.atoms[key]=Atom(key,candidate,tuple(query),reference_parts,parts,block_index,explained=explained)
        atom=self.atoms[key]
        if context:
            atom.contexts.add(context['context_id'])
            if key not in context['atom_ids']:context['atom_ids'].append(key)
        return atom

    def scan_once(self,index):
        if index in self.cache:return self.cache[index]
        b=self.chain[index];start=time.perf_counter();ids=[]
        for raw in b.parse_intra_alignment_sv(self.sequence,SVConfig.MIN_SV_SIZE,SVConfig.MERGE_INDEL_TOL,with_provenance=True):
            p=raw['evidence']['parts'];q=(min(x[1] for x in p),max(x[2] for x in p))
            c=SVCandidate(raw['type'],raw['chr'],raw['pos'],raw['len'],raw['read'],raw['info'])
            c.q_start,c.q_end=q
            a=self.atom(c,q,[[x[3],x[4]] for x in p],p,index)
            for context in self.block_contexts[index]:self.annotate_context(a,context)
            ids.append(a.atom_id)
        self.cigar_seconds+=time.perf_counter()-start
        self.cache[index]=ids;return ids

    def annotate_context(self,atom,context):
        atom.contexts.add(context['context_id'])
        if atom.atom_id not in context['atom_ids']:context['atom_ids'].append(atom.atom_id)

    def exact_overlap(self,a,b,start,end):
        if end<=start:return True
        aa=list(mapped_parts(a,(start,end),exact=True));bb=list(mapped_parts(b,(start,end),exact=True))
        if union_length((p[0],p[1]) for p in aa)!=end-start or union_length((p[0],p[1]) for p in bb)!=end-start:return False
        if end-start<self.k:return True
        return a.ref_name==b.ref_name and a.strand==b.strand and all(map_port(a,q)==map_port(b,q) for q in (start,end))

    def normalized_members(self,left,right,members,context):
        indices=([left] if left is not None else [])+members+([right] if right is not None else [])
        bounds={i:[self.chain[i].qs,self.chain[i].qe] for i in indices}
        for ia,ib in zip(indices,indices[1:]):
            a,b=self.chain[ia],self.chain[ib]
            overlap=max(0,bounds[ia][1]-bounds[ib][0])
            if overlap:
                start,end=bounds[ib][0],bounds[ia][1]
                if not self.exact_overlap(a,b,start,end):
                    context['issues'].append('QUERY_OWNERSHIP_CONFLICT');continue
                bounds[ib][0]=end
                context.setdefault('boundary_homology_bp',0);context['boundary_homology_bp']+=overlap
        if any(bounds[i][0]>=bounds[i][1] for i in members):context['issues'].append('EMPTY_BRANCH_AFTER_NORMALIZATION')
        return bounds

    def retained(self,chrom,start,end,exclude,extra=()):
        if end<=start:return False
        parts=self.source_index.get(chrom,[]);limit=bisect_left(self.source_starts.get(chrom,[]),end)
        begin=bisect_right(self.source_max_ends.get(chrom,[]),start,0,limit)
        ranges=[]
        for r0,r1,q0,q1,index in parts[begin:limit]:
            if r1<=start:continue
            # Clip the excluded query occurrence, rather than discarding a
            # complete long run because its boundary overlaps by a few bases.
            allowed=[(q0,min(q1,exclude[0])),(max(q0,exclude[1]),q1)]
            for qs,qe in allowed:
                if qe<=qs:continue
                if self.chain[index].strand=='+':a,b=r0+qs-q0,r0+qe-q0
                else:a,b=r1-(qe-q0),r1-(qs-q0)
                if a<end and b>start:ranges.append((max(start,a),min(end,b)))
        for qs,qe,rs,re,strand in extra:
            for a,b in ((qs,min(qe,exclude[0])),(max(qs,exclude[1]),qe)):
                if b<=a:continue
                lo,hi=(rs+a-qs,rs+b-qs) if strand=='+' else (re-(b-qs),re-(a-qs))
                if lo<end and hi>start:ranges.append((max(start,lo),min(end,hi)))
        return union_length(ranges)==end-start

    def emit_gap(self,q0,q1,chrom,r0,r1,strand,context,label='COMPLEX',explained=True):
        key=(q0,q1,chrom,r0,r1,strand)
        if key in self.gap_cache:
            answer=self.gap_cache[key]
            if answer:self.annotate_context(answer,context)
            return answer
        answer,handled=self.exact_linear_gap(q0,q1,chrom,r0,r1,strand,context)
        if not handled:
            answer,handled=self.padded_linear_gap(q0,q1,chrom,r0,r1,strand,context)
        if handled:
            self.gap_cache[key]=answer
            return answer
        if (1 if strand=='+' else -1)*(r1-r0)<=-SVConfig.MIN_SV_SIZE:
            context['issues'].append('INTRA_FRAME_DUPLICATION_UNRESOLVED')
            self.gap_cache[key]=None;return None
        if strand=='+':
            a=SimpleNamespace(qs=q0-1,qe=q0,ref_name=chrom,ref_start=r0-1,ref_end=r0,strand=strand)
            b=SimpleNamespace(qs=q1,qe=q1+1,ref_name=chrom,ref_start=r1,ref_end=r1+1,strand=strand)
        else:
            a=SimpleNamespace(qs=q0-1,qe=q0,ref_name=chrom,ref_start=r0,ref_end=r0+1,strand=strand)
            b=SimpleNamespace(qs=q1,qe=q1+1,ref_name=chrom,ref_start=r1-1,ref_end=r1,strand=strand)
        c=self.kernel(a,b,chrom,strand,label,isbackbone=False)
        answer=None
        if c:
            answer=self.atom(c,(min(q0,q1),max(q0,q1)),[[min(r0,r1),max(r0,r1)]],
                [['GAP',q0,q1,r0,r1,strand]],context=context,explained=explained)
        self.gap_cache[key]=answer;return answer

    def exact_linear_gap(self,q0,q1,chrom,r0,r1,strand,context):
        """Resolve the complete residual before using a net-length heuristic."""
        sigma=1 if strand=='+' else -1
        if q1<q0 or sigma*(r1-r0)<0 or self.reference is None:return None,False
        f=self.sequence[q0:q1].upper()
        r=self.refseq(chrom,*sorted((r0,r1)),strand)
        if r is None:return None,False
        if any(c not in 'ACGT' for c in f):
            context['issues'].append('ASSEMBLY_GAP');return None,True
        if any(c not in 'ACGT' for c in r):
            context['issues'].append('REFERENCE_GAP');return None,True
        model=exact_flank_model(f,r,k=self.k)
        if model is None:return None,False
        if not model.changes:return None,True
        if len(model.changes)!=1 or model.changes[0]['type']!='GAP':return None,False
        change=model.changes[0];qs,qe=(q0+x for x in change['query'])
        a,b=(r0+sigma*x for x in change['reference']);rs,re=sorted((a,b))
        size=max(qe-qs,re-rs)
        context.setdefault('sequence_checks',[]).append(dict(status='EXPLAINED',method='EXACT_FLANKS',
            query=[q0,q1],reference=[chrom,*sorted((r0,r1)),strand],homologies=model.homologies,derived_cigar=[]))
        if size<SVConfig.MIN_SV_SIZE:return None,True
        kind='INS' if qe>qs else 'DEL';info=dict(SOURCE='EXACT_GAP',RECHECK_METHOD='EXACT_FLANKS')
        if kind=='INS':info['SEQ']=f[change['query'][0]:change['query'][1]] if sigma==1 else reverse_complement(f[change['query'][0]:change['query'][1]])
        candidate=SVCandidate(kind,chrom,rs,size,self.name,info,end=re)
        candidate.q_start,candidate.q_end=qs,qe;candidate.reference_gap=(rs,re)
        return self.atom(candidate,(qs,qe),[[rs,re]],[['EXACT_GAP',qs,qe,rs,re,strand]],context=context),True

    def padded_linear_gap(self,q0,q1,chrom,r0,r1,strand,context):
        left,right=context.get('left_backbone'),context.get('right_backbone')
        if left is None or right is None:return None,False
        a,b=self.chain[left],self.chain[right];sigma=1 if strand=='+' else -1
        if a.ref_name!=chrom or b.ref_name!=chrom or a.strand!=strand or b.strand!=strand:return None,False
        pad=self.k-1+max(0,q0-q1,sigma*(r0-r1))
        if (q0-pad<a.qs or q1+pad>b.qe or
                union_length((s,e) for s,e,_,_ in mapped_parts(a,(q0-pad,q0),exact=True))!=pad or
                union_length((s,e) for s,e,_,_ in mapped_parts(b,(q1,q1+pad),exact=True))!=pad or
                map_port(a,q0-pad)!=r0-sigma*pad or map_port(b,q1+pad)!=r1+sigma*pad):
            return None,False
        atom,handled=self.exact_linear_gap(q0-pad,q1+pad,chrom,r0-sigma*pad,r1+sigma*pad,strand,context)
        if handled:
            context['gap_normalization']=dict(method='EXACT_PADDED_FLANKS',padding=pad,
                original_query=[q0,q1],original_reference=[chrom,r0,r1,strand])
        return atom,handled

    def sequence_gap(self,q0,q1,chrom,u0,u1,sigma,context):
        """Consume one query-oriented residual. Return atoms and completeness."""
        if q1<q0 or u1<u0:
            if max(q0-q1,u0-u1)<SVConfig.MIN_SV_SIZE:return [],True
            context['issues'].append('BOUNDARY_OVERLAP');return [],False
        strand='+' if sigma==1 else '-';r0,r1=sigma*u0,sigma*u1
        gq,gr=q1-q0,u1-u0
        if not gq and not gr:return [],True
        if min(gq,gr)<SVConfig.GAP_TOLERANCE:
            atom=self.emit_gap(q0,q1,chrom,r0,r1,strand,context,
                explained=min(gq,gr)<SVConfig.MIN_SV_SIZE)
            return ([atom] if atom else []),min(gq,gr)<SVConfig.MIN_SV_SIZE
        f=self.sequence[q0:q1].upper();r=self.refseq(chrom,min(r0,r1),max(r0,r1),strand)
        if r is None:context['issues'].append('REFERENCE_UNAVAILABLE');return [],False
        if context['frame']=='BACKBONE_GAP' and self.recover_inverse(f,r,q0,chrom,u0,sigma,context):return [],True
        result=check_sequence(f,r,k=self.k,max_bases=SVConfig.SEED_MAX_BASES,
            max_blocks=SVConfig.SEED_MAX_BLOCKS,min_gap=SVConfig.GAP_TOLERANCE,
            max_ref_gap=SVConfig.MIN_TRA_SIZE,realign=SVConfig.REPLACEMENT_REALIGN,preset=SVConfig.CALLER_PRESET)
        if context['frame']=='BACKBONE_GAP' and not result.issues and not result.homologies and len(result.changes)==1 and result.changes[0]['type']=='DELINS':
            recovered=self.padded_separator(q0,q1,chrom,u0,u1,sigma,context)
            if recovered is not None:return recovered,True
        context.setdefault('sequence_checks',[]).append(dict(check_id=identifier('CHECK',self.name,q0,q1,chrom,r0,r1),
            query=[q0,q1],reference=[chrom,min(r0,r1),max(r0,r1),strand],status=result.status,method=result.method,
            homologies=result.homologies,alternative_homologies=result.alternative_homologies,derived_cigar=result.alignment))
        context['issues'].extend(result.issues);atoms=[]
        for change in result.changes:
            qs,qe=(q0+v for v in change['query']);a,b=(u0+v for v in change['reference'])
            rs,re=sorted((sigma*a,sigma*b));kind=change['type']
            if kind=='GAP':
                atom=self.emit_gap(qs,qe,chrom,sigma*a,sigma*b,strand,context)
                if atom:atoms.append(atom)
                continue
            info=dict(SOURCE='SEQUENCE_CHECK',RECHECK_METHOD=result.method,RESOLUTION='RESOLVED',
                DECOMPOSITION='UNDETERMINED' if kind=='DELINS' else 'SINGLE')
            if kind in ('INS','DEL'):
                info['SOURCE']='ANCHOR_RECHECK';info['SEQ']=change.get('sequence','') if sigma==1 else reverse_complement(change.get('sequence',''))
                length=change['length']
            elif kind=='DELINS':
                info.update(INS_LEN=qe-qs,DEL_LEN=re-rs,NET_LEN=(qe-qs)-(re-rs),
                    SEQ=change['sequence'] if sigma==1 else reverse_complement(change['sequence']))
                length=abs((qe-qs)-(re-rs))
            else:length=re-rs
            if kind=='INV' and length<SVConfig.min_inv_size():continue
            c=SVCandidate(kind,chrom,rs,length,self.name,info,end=re)
            c.q_start,c.q_end=qs,qe
            atoms.append(self.atom(c,(qs,qe),[[rs,re]],[[kind,qs,qe,rs,re,strand]],context=context))
        return atoms,not result.issues

    def recover_inverse(self,f,r,q0,chrom,u0,sigma,context):
        model=exact_inverse_edit(f,r,k=self.k,min_edit=SVConfig.MIN_SV_SIZE,
            max_bases=SVConfig.EXACT_INVERSE_MAX_BASES,min_coverage=SVConfig.EXACT_INVERSE_MIN_COVERAGE)
        if model is None:return False
        source_start,source_end=sorted(sigma*(u0+x) for x in model['reference'])
        qspan=tuple(q0+x for x in model['query'])
        if any(self.inside(a,qspan,(source_start,source_end)) for a in self.atoms.values()):return False
        segments=[]
        for qs,qe,rs,re in model['inverse_arms']:
            lo,hi=sorted((sigma*(u0+rs),sigma*(u0+re)))
            segments.append(dict(block_index=None,query=[q0+qs,q0+qe],source=[chrom,lo,hi,'-' if sigma==1 else '+']))
        for i in context['member_block_indices']:
            for seg in segments:
                for s,e,a,b in mapped_parts(self.chain[i],seg['query']):
                    qs,qe=seg['query'];lo,hi=seg['source'][1:3]
                    expected=sorted((hi-(s-qs),hi-(e-qs))) if sigma==1 else (lo+s-qs,lo+e-qs)
                    if (a,b)!=tuple(expected):return False
        qs,qe=(q0+x for x in model['query_edit']);rs,re=sorted(sigma*(u0+x) for x in model['reference_edit'])
        ins,deleted=model['ins_len'],model['del_len'];kind='DELINS' if ins and deleted else 'INS' if ins else 'DEL'
        sequence=self.sequence[qs:qe]
        if sigma==1:sequence=reverse_complement(sequence)
        info=dict(SOURCE='EXACT_INVERSE',RECHECK_METHOD='EXACT_INVERSE',RESOLUTION='RESOLVED')
        if kind in ('INS','DELINS'):info['SEQ']=sequence
        if kind=='DELINS':info.update(INS_LEN=ins,DEL_LEN=deleted,NET_LEN=ins-deleted)
        c=SVCandidate(kind,chrom,rs,abs(ins-deleted),self.name,info,end=re);c.q_start,c.q_end=qs,qe
        atom=self.atom(c,(qs,qe),[[rs,re]],[[kind,qs,qe,rs,re,'-' if sigma==1 else '+']],context=context)
        context.update(frame='INV',source_intervals=[[chrom,source_start,source_end,'-' if sigma==1 else '+']],
            segments=segments,query_span=list(qspan),target_query_span=list(qspan),exact_inverse_model=model)
        context['issues']=[x for x in context['issues'] if x!='INVERSION_BOUNDARY_SOURCE_UNRESOLVED']
        self.frame_junctions(context,segments,chrom,sigma*(u0+model['reference'][0]),sigma*(u0+model['reference'][1]),'+' if sigma==1 else '-')
        self.make_parent('CPX',chrom,source_start,source_end,qspan,context,[atom],
            dict(FRAME='INV',INS_LEN=ins,DEL_LEN=deleted,INV_SPAN=source_end-source_start,NET_LEN=ins-deleted,RECHECK_METHOD='EXACT_INVERSE'))
        return True

    def padded_separator(self,q0,q1,chrom,u0,u1,sigma,context):
        pad=self.k-1;left,right=context['left_backbone'],context['right_backbone']
        if left is None or right is None:return None
        a,b=self.chain[left],self.chain[right]
        if q0-pad<a.qs or q1+pad>b.qe:return None
        if union_length((s,e) for s,e,_,_ in mapped_parts(a,(q0-pad,q0),exact=True))!=pad:return None
        if union_length((s,e) for s,e,_,_ in mapped_parts(b,(q1,q1+pad),exact=True))!=pad:return None
        if map_port(a,q0-pad)!=sigma*(u0-pad) or map_port(b,q1+pad)!=sigma*(u1+pad):return None
        f=self.sequence[q0-pad:q1+pad].upper();r=self.refseq(chrom,*sorted((sigma*(u0-pad),sigma*(u1+pad))),'+' if sigma==1 else '-')
        if r is None:return None
        check=check_sequence(f,r,k=self.k,max_bases=SVConfig.SEED_MAX_BASES,max_blocks=SVConfig.SEED_MAX_BLOCKS,
            min_gap=SVConfig.GAP_TOLERANCE,max_ref_gap=SVConfig.MIN_TRA_SIZE)
        if check.issues or len(check.changes)!=2 or any(c['type']!='GAP' for c in check.changes):return None
        trimmed=[]
        for change in check.changes:
            qs,qe=change['query'];rs,re=change['reference'];x=common_prefix(f[qs:qe],r[rs:re]);qs+=x;rs+=x
            y=common_prefix(f[qs:qe][::-1],r[rs:re][::-1]);qe-=y;re-=y
            if min(qe-qs,re-rs)!=0 or max(qe-qs,re-rs)<SVConfig.MIN_SV_SIZE+SVConfig.INSDEL_TOLERANCE:return None
            trimmed.append((qs,qe,rs,re))
        if sorted('INS' if qe>qs else 'DEL' for qs,qe,rs,re in trimmed)!=['DEL','INS']:return None
        # An aligner may absorb the last island base into the right anchor as
        # well as the first base into the left anchor. The separator must
        # intersect both original gaps; the exact flanking proof and the two
        # SV-sized pure residuals certify the padded representation.
        if not any(e-s>=self.k and e>pad and s<len(f)-pad and re>pad and rs<len(r)-pad
                   for s,e,rs,re in check.homologies):return None
        atoms=[]
        for qs,qe,rs,re in trimmed:
            atom=self.emit_gap(q0-pad+qs,q0-pad+qe,chrom,sigma*(u0-pad+rs),sigma*(u0-pad+re),'+' if sigma==1 else '-',context)
            if atom:atom.candidate.info['RECHECK_METHOD']='EXACT_SEPARATOR';atoms.append(atom)
        context.setdefault('sequence_checks',[]).append(dict(status='EXPLAINED',method='EXACT_SEPARATOR',
            query=[q0-pad,q1+pad],reference=[chrom,*sorted((sigma*(u0-pad),sigma*(u1+pad)))],homologies=check.homologies,derived_cigar=[]))
        context['separators'].extend(dict(query=[q0-pad+s,q0-pad+e],reference=[chrom,*sorted((sigma*(u0-pad+rs),sigma*(u0-pad+re)))],origin='EXACT_PADDED_CHECK') for s,e,rs,re in check.homologies)
        return atoms

    def ordinary(self,left,right):
        a,b=self.chain[left],self.chain[right];context=self.context(left,right,[],'BACKBONE_GAP')
        if a.ref_name==b.ref_name and a.strand==b.strand:
            direction=1 if a.strand=='+' else -1
            if direction*(entry(b)-exit_port(a))<=-SVConfig.MIN_SV_SIZE:
                self.tandem_backbone(context,left,right)
                self.scan_once(left);self.scan_once(right)
                return
            if self.exact_backbone_inversion(context,left,right):
                self.scan_once(left);self.scan_once(right)
                return
            if self.exact_backbone_indel(context,left,right):
                self.scan_once(left);self.scan_once(right)
                return
        # Use the actual original blocks for the original tandem overlap guard.
        c=self.kernel(a,b,a.ref_name,a.strand,'backbone',isbackbone=False)
        if c:
            atom=self.atom(c,(min(a.qe,b.qs),max(a.qe,b.qs)),
                [[min(exit_port(a),entry(b)),max(exit_port(a),entry(b))]],
                [['BACKBONE',left,right,a.qe,b.qs,exit_port(a),entry(b)]],context=context)
            context.update(resolution='RESOLVED',decomposition='SINGLE',frame='LOCAL_INDEL')
        elif a.ref_name==b.ref_name and a.strand==b.strand:
            sigma=1 if a.strand=='+' else -1
            if b.qs-a.qe>=SVConfig.GAP_TOLERANCE and sigma*(entry(b)-exit_port(a))>=SVConfig.GAP_TOLERANCE:
                _,complete=self.sequence_gap(a.qe,b.qs,a.ref_name,sigma*exit_port(a),sigma*entry(b),sigma,context)
                if not context.get('event_id'):context.update(resolution='RESOLVED' if complete else 'PARTIAL',decomposition='SEPARABLE' if len(context['atom_ids'])>1 else 'UNDETERMINED')
            else:context.update(resolution='RESOLVED',decomposition='SINGLE',frame='COLLINEAR_NO_CALL')
        else:context['issues'].append('INCOMPATIBLE_BACKBONE')
        self.scan_once(left);self.scan_once(right)

    def exact_backbone_indel(self,context,left,right):
        a,b=self.chain[left],self.chain[right]
        atom,handled=self.padded_linear_gap(a.qe,b.qs,a.ref_name,exit_port(a),entry(b),a.strand,context)
        if not handled:atom,handled=self.exact_linear_gap(a.qe,b.qs,a.ref_name,exit_port(a),entry(b),a.strand,context)
        if not handled:return False
        context.update(resolution='UNRESOLVED' if context['issues'] else 'RESOLVED',
            decomposition='SINGLE',frame='LOCAL_INDEL' if atom else 'COLLINEAR_NO_CALL')
        return True

    def exact_backbone_inversion(self,context,left,right):
        """Check short and long inverse residuals before net-gap dispatch.

        Exact padding is taken only from owned mapped anchor bases. This can
        recover an equivalent inversion boundary absorbed into either anchor.
        """
        a,b=self.chain[left],self.chain[right];sigma=1 if a.strand=='+' else -1
        q0,q1=a.qe,b.qs;u0,u1=sigma*exit_port(a),sigma*entry(b)
        if q1<=q0 or u1<=u0:return False
        pad=self.k-1
        if (q0-pad<a.qs or q1+pad>b.qe or
                union_length((s,e) for s,e,_,_ in mapped_parts(a,(q0-pad,q0),exact=True))!=pad or
                union_length((s,e) for s,e,_,_ in mapped_parts(b,(q1,q1+pad),exact=True))!=pad or
                map_port(a,q0-pad)!=sigma*(u0-pad) or map_port(b,q1+pad)!=sigma*(u1+pad)):
            pad=0
        f=self.sequence[q0-pad:q1+pad].upper()
        r=self.refseq(a.ref_name,*sorted((sigma*(u0-pad),sigma*(u1+pad))),a.strand)
        if r is None:return False
        model=exact_flank_model(f,r,k=self.k)
        if model is None or len(model.changes)!=1 or model.changes[0]['type']!='INV':return False
        change=model.changes[0];s,t=change['reference']
        rs,re=sorted((sigma*(u0-pad+s),sigma*(u0-pad+t)))
        qs,qe=(q0-pad+x for x in change['query'])
        source=[a.ref_name,rs,re,'-' if sigma==1 else '+']
        segment=dict(block_index=None,query=[qs,qe],source=source)
        context.update(frame='INV',target_strand=a.strand,source_intervals=[source],
            segments=[segment],query_span=[qs,qe],target_query_span=[qs,qe])
        context.setdefault('sequence_checks',[]).append(dict(status='EXPLAINED',method='EXACT_FLANKS',
            query=[q0-pad,q1+pad],reference=[a.ref_name,*sorted((sigma*(u0-pad),sigma*(u1+pad))),a.strand],
            homologies=model.homologies,derived_cigar=[],padding=pad))
        self.frame_junctions(context,[segment],a.ref_name,sigma*(u0-pad+s),sigma*(u0-pad+t),a.strand)
        self.make_parent('INV',a.ref_name,rs,re,(qs,qe),context,
            info=dict(FRAME='INV',NET_LEN=0,INV_SPAN=re-rs,RECHECK_METHOD='EXACT_FLANKS'))
        return True

    def claim(self,atom,parent,context,copy=False):
        if atom.parent and atom.parent!=parent.id:
            context['issues'].append('ATOM_OWNERSHIP_CONFLICT');return False
        atom.parent=parent.id;self.annotate_context(atom,context)
        atom.candidate.info.update(ROLE='COMPONENT',PARENT_EVENT=parent.id,
            ALLELE_CONTEXT='ADDED_DUPLICATE' if copy else 'REFERENCE')
        return True

    def insertion_signature(self,atom):
        if atom.candidate.type not in ('INS','DELINS'):return ''
        if atom.block_index is not None:
            sequence=''.join(self.sequence[p[1]:p[2]] if self.chain[atom.block_index].strand=='+' else
                reverse_complement(self.sequence[p[1]:p[2]]) for p in atom.parts if p[0]==1)
        else:sequence=atom.candidate.info.get('SEQ','')
        return hashlib.sha256(sequence.upper().encode()).hexdigest()

    def make_parent(self,kind,chrom,start,end,qspan,context,children=(),info=None,basis=None,basis_reference=False):
        if kind=='INV' and not children and len(context.get('source_intervals',[]))==1:
            from asmsv_next.allele_normalization import canonical_inversion
            orientation=context.get('target_strand',self.strand)
            normalized=canonical_inversion(self.reference,chrom,start,end,self.sequence,qspan,orientation,self.k-1)
            if normalized is not None:
                old_start,old_end,old_query=start,end,qspan
                start,end,qspan=normalized
                source=[chrom,start,end,'-' if orientation=='+' else '+']
                context['original_source_intervals']=context['source_intervals']
                context['source_intervals']=[source]
                context['original_junctions']=context.get('junctions',[])
                context['query_span']=context['target_query_span']=list(qspan)
                context['inversion_normalization']=dict(method='EXACT_BOUNDARY_EQUIVALENCE',
                    original_reference=[old_start,old_end],original_query=list(old_query),
                    reference=[start,end],query=list(qspan))
                sigma=1 if orientation=='+' else -1
                first,last=(start,end) if sigma==1 else (end,start)
                self.frame_junctions(context,[dict(block_index=None,query=list(qspan),source=source)],chrom,first,last,orientation)
                for edge in context['junctions']:edge['evidence']='ALLELE_EQUIVALENT_PORTS'
                info=dict(info or {},INV_SPAN=end-start)
        if kind=='INV' and end-start<SVConfig.min_inv_size():
            context.update(frame='BELOW_MIN_INV',resolution='RESOLVED',decomposition='NONE',junctions=[],
                size_policy=dict(kind='INV',length=end-start,minimum=SVConfig.min_inv_size()))
            return None
        metadata=dict(SOURCE='STRUCTURE',RESOLUTION='RESOLVED',DECOMPOSITION='COUPLED' if kind=='CPX' else 'SINGLE',
            CONTEXT_ID=context['context_id'],ROLE='EVENT',EVIDENCE_SCOPE='SELECTED_PATH',ALT_ASSESSMENT='NOT_AVAILABLE')
        metadata.update(info or {})
        frame=metadata.get('FRAME');sources=[list(s) for s in context.get('source_intervals',[])]
        if frame in ('INV','DUP'):
            for source in sources:source[3]='-' if frame=='INV' or metadata.get('DUP_ORIENT')=='INVERTED' else '+'
        elif frame=='MULTISOURCE' and context.get('target_strand')=='-':
            sources=[s[:3]+['-' if s[3]=='+' else '+'] for s in reversed(sources)]
            for key in ('CPX_COMPONENTS','CPX_INTERVALS','SOURCE_STATES'):
                if metadata.get(key):metadata[key]=','.join(reversed(metadata[key].split(',')))
        if frame in ('DUP','MULTISOURCE'):metadata.setdefault('SOURCE_SEARCH_SCOPE','SELECTED_PATH')
        metadata['SOURCE_SIGNATURE']='|'.join(':'.join(map(str,s)) for s in sources)
        if kind=='CPX':
            metadata.pop('CSV_TYPE',None);metadata.pop('SUBTYPE',None)
            if metadata.get('FRAME')=='MULTISOURCE':
                metadata.setdefault('CPX_TYPE','COMPLEX')
            else:metadata.update(complex_metadata(frame,sources,children))
        metadata['STRUCTURE_KEY']=identifier('STRUCT',sources,
            [metadata.get('TARGET_CHR'),metadata.get('TARGET_POS'),metadata.get('DUP_ORIENT'),metadata.get('DUP_MODE')],
            sorted((a.candidate.type,a.candidate.chrom,a.candidate.pos,a.candidate.end,a.candidate.sv_len,self.insertion_signature(a)) for a in children))
        event_id=identifier('EVENT',self.haplotype,self.name,context['context_id'],kind,start,end,metadata.get('CPX_TYPE'),metadata.get('CPX_COMPONENTS'))
        metadata['EVENT']=event_id
        size=end-start if kind=='INV' else abs(int(metadata.get('NET_LEN',end-start)))
        parent=SVCandidate(kind,chrom,start,size,self.name,metadata,end=end)
        if frame=='DUP' and SVConfig.MAX_SV_SIZE!=-1 and metadata.get('NET_LEN',0)>=SVConfig.MAX_SV_SIZE:
            parent.filter='OVERSIZE'
        parent.id=event_id;parent.q_start,parent.q_end=qspan
        if any(atom.parent and atom.parent!=event_id for atom in list(children)+([basis] if basis else [])):
            context['issues'].append('ATOM_OWNERSHIP_CONFLICT');return None
        for atom in children:
            if not self.claim(atom,parent,context,metadata.get('FRAME')=='DUP'):return None
        if basis is not None:self.claim(basis,parent,context,metadata.get('FRAME')=='DUP' and not basis_reference)
        connections=context.setdefault('junctions',[])
        for atom in children:
            for part in atom.parts:
                if part[0] not in (1,2):continue
                op,qs,qe,rs,re=part
                left=dict(chrom=atom.candidate.chrom,position=rs,side='R')
                right=dict(chrom=atom.candidate.chrom,position=re,side='L')
                edge=junction(self.haplotype,self.name,context['context_id'],left,right,[qs,qe],atom.atom_id if op==1 else None)
                if all(x['junction_id']!=edge['junction_id'] for x in connections):connections.append(edge)
        self.parents.append(parent)
        context.update(resolution='RESOLVED',decomposition=metadata['DECOMPOSITION'],event_id=event_id)
        self.events.append(dict(event_id=event_id,svtype=kind,cpx_type=metadata.get('CPX_TYPE'),
            context_id=context['context_id'],query_name=self.name,haplotype=self.haplotype,
            target=[metadata.get('TARGET_CHR',chrom),metadata.get('TARGET_POS',start),metadata.get('TARGET_POS',end)],
            reference_span=[chrom,start,end],query_span=list(qspan),components=[a.atom_id for a in children],
            basis_atom=basis.atom_id if basis else None,info={k:v for k,v in metadata.items() if k!='SEQ'},
            segments=context.get('segments',[]),junctions=context.get('junctions',[])))
        return parent

    def insertion_projection(self,parent,chrom,target,qspan,context,strand):
        """Reference-change view of a DUP; never an additional native event."""
        qs,qe=qspan
        if qe<=qs:return
        sequence=self.sequence[qs:qe]
        if strand=='-':sequence=reverse_complement(sequence)
        candidate=SVCandidate('INS',chrom,target,qe-qs,self.name,
            dict(SOURCE='DUP_PROJECTION',ROLE='PROJECTION',PARENT_EVENT=parent.id,CONTEXT_ID=context['context_id'],
                 SEQ=sequence,EVIDENCE_SCOPE='SELECTED_PATH'),filter_str=parent.filter)
        candidate.q_start,candidate.q_end=qs,qe
        candidate.id=identifier('PROJ',parent.id,'INS',chrom,target,qs,qe)
        self.projections.append(candidate)
        context.setdefault('insertion_projections',[]).append(dict(projection_id=candidate.id,parent_event=parent.id,
            target=[chrom,target],query=[qs,qe],length=qe-qs,sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest()))

    def tandem_backbone(self,context,left,right):
        """Two reference visits, independently supported on the physical query.

        Size never selects INS versus DUP. A positive query filler is first
        tested as an exact extension; otherwise it is explicit junction sequence.
        """
        a,b=self.chain[left],self.chain[right];sigma=1 if a.strand=='+' else -1
        u0,u1=sigma*entry(b),sigma*exit_port(a)
        context.update(frame='DUP',target_strand=a.strand)
        if not sigma*entry(a)<=u0<u1<=sigma*exit_port(b):
            context['issues'].append('NON_TANDEM_REFERENCE_ORDER');return
        qa,qb=a.qe,b.qs;gap=qb-qa;extra=[];normalization='NONE'
        if gap<0:
            if -gap>=self.k or not self.exact_overlap(a,b,qb,qa):
                context['issues'].append('QUERY_OWNERSHIP_CONFLICT');return
            point=map_port(b,qa)
            if point is None:context['issues'].append('UNMAPPED_DUP_PORT');return
            u0=sigma*point;normalization='MICROHOMOLOGY';qb=qa
        elif gap>0:
            if gap>MAX_LINEAR_BASES:
                context['issues'].append('DUP_JUNCTION_BUDGET');return
            filler=self.sequence[qa:qb].upper()
            if any(c not in 'ACGT' for c in filler):
                context['issues'].append('ASSEMBLY_GAP');return
            rs,re=sorted((sigma*(u0-gap),sigma*u0))
            if sigma*entry(a)<=u0-gap and filler==self.refseq(a.ref_name,rs,re,a.strand):
                extra.append((qa,qb,rs,re,a.strand));u0-=gap;normalization='RIGHT_ANCHOR_EXTENSION';qb=qa
            else:
                rs,re=sorted((sigma*u1,sigma*(u1+gap)))
                if u1+gap<=sigma*exit_port(b) and filler==self.refseq(a.ref_name,rs,re,a.strand):
                    extra.append((qa,qb,rs,re,a.strand));u1+=gap;normalization='LEFT_ANCHOR_EXTENSION';qa=qb
        s,t=sorted((sigma*u0,sigma*u1));length=t-s
        context['normalization']=normalization
        if extra:context['boundary_extensions']=[dict(query=[q0,q1],reference=[a.ref_name,r0,r1,orientation],method='EXACT')
            for q0,q1,r0,r1,orientation in extra]
        if length<SVConfig.MIN_SV_SIZE:
            context.update(frame='BELOW_MIN_DUP',resolution='RESOLVED',decomposition='SINGLE');return
        occurrences=[]
        for index in (left,right):
            mapped=list(mapped_parts(self.chain[index],reference=(s,t)))
            if not mapped:context['issues'].append('DUP_SOURCE_UNMAPPED');return
            span=[min(p[0] for p in mapped),max(p[1] for p in mapped)]
            if extra and ((index==right and normalization=='RIGHT_ANCHOR_EXTENSION') or
                          (index==left and normalization=='LEFT_ANCHOR_EXTENSION')):
                span=[min(span[0],extra[0][0]),max(span[1],extra[0][1])]
            if index==right and normalization=='MICROHOMOLOGY':span[0]=max(span[0],qb)
            occurrences.append(span)
        if occurrences[0][1]>occurrences[1][0]:
            context['issues'].append('QUERY_OWNERSHIP_CONFLICT');return
        if min(e-s for s,e in occurrences)<SVConfig.MIN_SV_SIZE:
            context['issues'].append('DUP_QUERY_SUPPORT_TOO_SHORT');return
        context.update(source_intervals=[[a.ref_name,s,t,a.strand]],
            query_span=[occurrences[0][0],occurrences[1][1]],
            segments=[dict(block_index=i,query=list(span),source=[a.ref_name,s,t,a.strand])
                      for i,span in zip((left,right),occurrences)])
        # A framework exists before scanning its original CIGAR views.
        for i in (left,right):self.scan_once(i)
        edits=[[self.atoms[k] for k in self.cache[index] if self.inside(self.atoms[k],span,(s,t))]
               for index,span in zip((left,right),occurrences)]
        if edits[0] and edits[1]:context['issues'].append('DUP_OCCURRENCE_UNRESOLVED');return
        affected=0 if edits[0] else 1 if edits[1] else (1 if a.strand=='+' else 0)
        qs,qe=occurrences[affected];children=edits[affected]
        if not self.retained(a.ref_name,s,t,(qs,qe),extra):
            context['issues'].append('DUP_SOURCE_RETENTION_UNRESOLVED');return
        if not self.safe_children(children,context):return
        # An edited left genomic occurrence is inserted before the unedited one.
        target=s if (affected==0)==(a.strand=='+') else t
        jq0,jq1=occurrences[0][1],occurrences[1][0];junction_length=jq1-jq0
        added=(occurrences[0][0],occurrences[1][0]) if affected==0 else (occurrences[0][1],occurrences[1][1])
        if any(c not in 'ACGT' for c in self.sequence[added[0]:added[1]].upper()):
            context['issues'].append('ASSEMBLY_GAP');return
        if junction_length>=SVConfig.MIN_SV_SIZE:
            sequence=self.sequence[jq0:jq1]
            if a.strand=='-':sequence=reverse_complement(sequence)
            c=SVCandidate('INS',a.ref_name,target,junction_length,self.name,
                dict(SOURCE='DUP_JUNCTION',SEQ=sequence,ALLELE_CONTEXT='DUPLICATION_JUNCTION'))
            c.q_start,c.q_end=jq0,jq1
            children=children+[self.atom(c,(jq0,jq1),[[target,target]],
                [['DUP_JUNCTION',jq0,jq1,target,target,a.strand]],context=context)]
        info=duplication_metadata(a.ref_name,s,t,a.ref_name,target,'DIRECT',precise=junction_length==0)
        info.update(FRAME='DUP',SOURCE_STATE='RETAINED',DUP_TEMPLATE_LEN=length,
            DUP_SEQUENCE_LEN=qe-qs,DUP_OCCURRENCE=affected,NET_LEN=added[1]-added[0],
            JUNCTION_INS_LEN=junction_length,DUP_NORMALIZATION=normalization)
        ins,deleted=self.edit_lengths(children);info.update(INS_LEN=ins,DEL_LEN=deleted)
        p=dict(chrom=a.ref_name,position=sigma*u1,side='R' if sigma==1 else 'L')
        q=dict(chrom=a.ref_name,position=sigma*u0,side='L' if sigma==1 else 'R')
        context['junctions']=[junction(self.haplotype,self.name,context['context_id'],p,q,[jq0,jq1])]
        context['added_query_span']=list(added)
        kind='CPX' if children else 'DUP'
        parent=self.make_parent(kind,a.ref_name,s,t,(occurrences[0][0],occurrences[1][1]),context,children,info)
        if parent:
            for atom in children:
                if atom.candidate.info.get('SOURCE')=='DUP_JUNCTION':atom.candidate.info['ALLELE_CONTEXT']='DUPLICATION_JUNCTION'
            if SVConfig.MAX_SV_SIZE!=-1 and added[1]-added[0]>=SVConfig.MAX_SV_SIZE:parent.filter='OVERSIZE'
            self.insertion_projection(parent,a.ref_name,target,added,context,a.strand)

    @staticmethod
    def inside(atom,qspan,rspan):
        if not qspan[0]<=atom.query[0]<=atom.query[1]<=qspan[1]:return False
        return all(rspan[0]<=a<=b<=rspan[1] for a,b in atom.reference_parts)

    @staticmethod
    def edit_lengths(atoms):
        ins=sum(a.candidate.sv_len if a.candidate.type=='INS' else int(a.candidate.info.get('INS_LEN',0)) for a in atoms)
        deleted=sum(a.candidate.sv_len if a.candidate.type=='DEL' else int(a.candidate.info.get('DEL_LEN',0)) for a in atoms)
        return ins,deleted

    @staticmethod
    def safe_children(atoms,context):
        if any(a.candidate.type not in ('INS','DEL','DELINS') or not a.explained for a in atoms):
            context['issues'].append('COMPONENTS_UNEXPLAINED');return False
        return True

    def residual(self,q0,q1,chrom,u0,u1,sigma,side,context):
        """Restore a complete filler before the sole remaining-gap extraction."""
        if q1>=q0 and u1>=u0 and q1>q0:
            f=self.sequence[q0:q1].upper();r=self.refseq(chrom,min(sigma*u0,sigma*u1),max(sigma*u0,sigma*u1),'+' if sigma==1 else '-')
            if r is not None and f==r:return [],True,False
            if r is not None and exact_filler(f,r,side,self.k):
                if side=='right':
                    context['separators'].append(dict(query=[q0,q1],reference=[chrom,*sorted((sigma*u0,sigma*(u0+len(f))))],origin='EXACT_FILLER_CHECK'))
                    atoms,ok=self.sequence_gap(q1,q1,chrom,u0+len(f),u1,sigma,context)
                else:
                    context['separators'].append(dict(query=[q0,q1],reference=[chrom,*sorted((sigma*(u1-len(f)),sigma*u1))],origin='EXACT_FILLER_CHECK'))
                    atoms,ok=self.sequence_gap(q0,q0,chrom,u0,u1-len(f),sigma,context)
                return atoms,ok,True
        atoms,ok=self.sequence_gap(q0,q1,chrom,u0,u1,sigma,context)
        return atoms,ok,False

    def inversion_boundary_check(self,left_query,right_query,chrom,left_ref,right_ref,sigma,context):
        """An unaligned inverted source must not become a novel INS plus a DEL.

        A shared exact MIN_SV_SIZE-mer across opposite residuals is sufficient
        to withhold a complete frame. Repeats may therefore remain unresolved.
        This check changes grouping only; the original indel atoms are retained.
        """
        for query,own_ref,opposite_query,opposite_ref in (
                (left_query,left_ref,right_query,right_ref),
                (right_query,right_ref,left_query,left_ref)):
            nq=query[1]-query[0];nr=opposite_ref[1]-opposite_ref[0]
            if min(nq,nr)<SVConfig.MIN_SV_SIZE:continue
            if own_ref[1]-own_ref[0]>=SVConfig.MIN_SV_SIZE or opposite_query[1]-opposite_query[0]>=SVConfig.MIN_SV_SIZE:continue
            if nq+nr>SVConfig.SEED_MAX_BASES:
                context['issues'].append('INVERSION_BOUNDARY_RECHECK_BUDGET');continue
            f=self.sequence[query[0]:query[1]].upper()
            r=self.refseq(chrom,*sorted((sigma*opposite_ref[0],sigma*opposite_ref[1])),
                '+' if sigma==1 else '-')
            if r is None:
                context['issues'].append('REFERENCE_UNAVAILABLE');continue
            if shared_seed(f,reverse_complement(r),max(self.k,SVConfig.MIN_SV_SIZE)):
                context['issues'].append('INVERSION_BOUNDARY_SOURCE_UNRESOLVED')

    def branch(self,left,right,members):
        if left is not None and right is not None:
            a,b=self.chain[left],self.chain[right]
            if a.ref_name==b.ref_name and a.strand==b.strand:
                sigma=1 if a.strand=='+' else -1;lo=sigma*exit_port(a);hi=sigma*entry(b);cuts=[]
                for i in members:
                    h=self.chain[i]
                    if h.ref_name!=a.ref_name or h.strand!=a.strand:continue
                    if not lo<=sigma*entry(h)<sigma*exit_port(h)<=hi:continue
                    if not any(e-s>=self.k for s,e,_,_ in mapped_parts(h,exact=True)):continue
                    good=True
                    for j in members:
                        if i==j:continue
                        x=self.chain[j]
                        if x.ref_name!=a.ref_name or x.strand==a.strand:
                            good=False;break
                        xlo,xhi=sorted((sigma*x.ref_start,sigma*x.ref_end))
                        if (j<i and xhi>sigma*entry(h)) or (j>i and xlo<sigma*exit_port(h)):
                            good=False;break
                    if good:cuts.append(i)
                if cuts:
                    for x,y in zip([left]+cuts,cuts+[right]):
                        if y==x+1:self.ordinary(x,y)
                        else:self.branch(x,y,list(range(x+1,y)))
                    return
        context=self.context(left,right,members,'BRANCH')
        if left is None or right is None:
            context['frame']='OPEN'
            self.open_junction(context,left,right,members);return
        a,b=self.chain[left],self.chain[right];sigma=1 if a.strand=='+' else -1
        context['target_strand']=a.strand
        bounds=self.normalized_members(left,right,members,context)
        if context['issues']:return
        if a.ref_name!=b.ref_name or a.strand!=b.strand:
            context['issues'].append('INCOMPATIBLE_TARGET_ANCHORS');return
        branches=[self.chain[i] for i in members]
        descriptions=[]
        for i in members:
            block=self.chain[i];qs,qe=bounds[i]
            x,y=map_port(block,qs),map_port(block,qe)
            if x is None or y is None:context['issues'].append('UNMAPPED_BRANCH_PORT');return
            descriptions.append(dict(block_index=i,query=[qs,qe],source=[block.ref_name,min(x,y),max(x,y),block.strand]))
        context['segments']=descriptions
        q0,q1=bounds[left][1],bounds[right][0]
        context['target_query_span']=[q0,q1]
        port0,port1=map_port(a,bounds[left][1]),map_port(b,bounds[right][0])
        if port0 is None or port1 is None:
            context['issues'].append('UNMAPPED_TARGET_PORT');return
        target0,target1=sigma*port0,sigma*port1
        sources=[d['source'] for d in descriptions]
        same_target=all(s[0]==a.ref_name for s in sources)
        opposite=all(x.strand!=a.strand for x in branches)
        slo=min(min(sigma*s[1],sigma*s[2]) for s in sources) if same_target else None
        shi=max(max(sigma*s[1],sigma*s[2]) for s in sources) if same_target else None
        internal=[];complete=True
        inv_candidate=same_target and opposite and target0-self.k<slo<shi<target1+self.k
        # Geometry defines the frame before the sole original CIGAR SV scan.
        context['frame']='INV' if inv_candidate else 'SOURCE_PATH'
        for i in [left]+members+[right]:self.scan_once(i)
        for i,j in (zip(members,members[1:]) if inv_candidate else []):
            x,y=self.chain[i],self.chain[j]
            if x.ref_name==y.ref_name and x.strand==y.strand:
                direction=1 if x.strand=='+' else -1
                atoms,ok=self.sequence_gap(bounds[i][1],bounds[j][0],x.ref_name,
                    direction*map_port(x,bounds[i][1]),direction*map_port(y,bounds[j][0]),direction,context)
                internal.extend(atoms);complete=complete and ok
        if inv_candidate:
            context['frame']='INV'
            # All inverted source steps must move in the same reverse target order.
            if any(sigma*map_port(self.chain[i],bounds[i][1])<sigma*map_port(self.chain[j],bounds[j][0])-self.k for i,j in zip(members,members[1:])):
                context['issues'].append('SOURCE_ORDER_CONFLICT');return
            ls,le=sorted((sigma*slo,sigma*shi))
            self.inversion_boundary_check((q0,bounds[members[0]][0]),(bounds[members[-1]][1],q1),
                a.ref_name,(target0,slo),(shi,target1),sigma,context)
            if not internal and context['issues']==['INVERSION_BOUNDARY_SOURCE_UNRESOLVED']:
                f=self.sequence[q0:q1].upper();r=self.refseq(a.ref_name,*sorted((sigma*target0,sigma*target1)),a.strand)
                if r is not None and self.recover_inverse(f,r,q0,a.ref_name,target0,sigma,context):return
            left_atoms,lc,lsep=self.residual(q0,bounds[members[0]][0],a.ref_name,target0,slo,sigma,'left',context)
            right_atoms,rc,rsep=self.residual(bounds[members[-1]][1],q1,a.ref_name,shi,target1,sigma,'right',context)
            children=list(internal)
            for i in members:children.extend(self.atoms[k] for k in self.cache[i])
            if not lsep:children.extend(left_atoms)
            if not rsep:children.extend(right_atoms)
            children=list({x.atom_id:x for x in children}.values())
            complete=complete and lc and rc and not context['issues']
            context['source_intervals']=[[a.ref_name,ls,le,branches[0].strand]]
            if not complete or not self.safe_children(children,context):
                context['issues'].append('COMPONENTS_UNEXPLAINED')
                # The inversion itself remains observable even when the larger grouping is partial.
                if le-ls>=SVConfig.min_inv_size():
                    c=SVCandidate('INV',a.ref_name,ls,le-ls,self.name,dict(SOURCE='STRUCTURE',CONTEXT_ID=context['context_id'],RESOLUTION='PARTIAL'),filter_str='PARTIAL')
                    self.atom(c,(bounds[members[0]][0],bounds[members[-1]][1]),[[ls,le]],[['INV',ls,le]],context=context)
                return
            p0=slo if lsep or not left_atoms else target0;p1=shi if rsep or not right_atoms else target1
            start,end=sorted((sigma*p0,sigma*p1));ins,deleted=self.edit_lengths(children)
            kind='CPX' if children else 'INV'
            extra=dict(FRAME='INV',INS_LEN=ins,DEL_LEN=deleted,INV_SPAN=le-ls,NET_LEN=ins-deleted if children else 0)
            self.frame_junctions(context,descriptions,a.ref_name,sigma*p0,sigma*p1,a.strand)
            event_query=(min([bounds[members[0]][0]]+[x.query[0] for x in children]),
                         max([bounds[members[-1]][1]]+[x.query[1] for x in children]))
            self.make_parent(kind,a.ref_name,start,end,event_query,context,children,extra)
            if lsep or rsep:context['decomposition']='SEPARABLE'
            return
        # Source insertions: a closed target with no SV-sized target deletion.
        if abs(target1-target0)>=SVConfig.MIN_SV_SIZE:
            context['issues'].append('TARGET_REPLACEMENT_WITH_SOURCE_PATH');return
        if bounds[members[0]][0]-q0>=SVConfig.MIN_SV_SIZE or q1-bounds[members[-1]][1]>=SVConfig.MIN_SV_SIZE:
            context['issues'].append('INSERTION_FILLER_UNEXPLAINED');return
        runs=[]
        for d in descriptions:
            if runs:
                prev=runs[-1];s=d['source'];p=prev['source']
                direction=1 if s[3]=='+' else -1
                rg=direction*((s[1] if direction==1 else s[2])-(p[2] if direction==1 else p[1]))
                qg=d['query'][0]-prev['query'][1]
                if s[0]==p[0] and s[3]==p[3] and 0<=rg<SVConfig.MIN_SV_SIZE and 0<=qg<SVConfig.MIN_SV_SIZE:
                    prev['query'][1]=d['query'][1];p[1],p[2]=min(p[1],s[1]),max(p[2],s[2]);continue
            runs.append(dict(query=list(d['query']),source=list(d['source'])))
        context['source_intervals']=[d['source'] for d in runs]
        children=list(internal)
        for i in members:children.extend(self.atoms[k] for k in self.cache[i])
        children=list({x.atom_id:x for x in children}.values())
        if len(runs)==1:
            run=runs[0];chrom,s,t,orientation=run['source'];qspan=run['query']
            retained=self.retained(chrom,s,t,qspan)
            state='RETAINED' if retained else 'UNKNOWN'
            context['frame']='DUP' if retained else 'TEMPLATED_INSERTION'
            if children and (not retained or not self.safe_children(children,context) or not complete):
                context['issues'].append('DUP_EDIT_UNRESOLVED');return
            direct='DIRECT' if orientation==a.strand else 'INVERTED'
            target=min(sigma*target0,sigma*target1)
            filler=(bounds[members[0]][0]-q0)+(q1-bounds[members[-1]][1])
            info=dict(SOURCE_STATE=state,FRAME='DUP' if retained else 'INS',SOURCE_CHR=chrom,SOURCE_START=s,SOURCE_END=t,
                TARGET_CHR=a.ref_name,TARGET_POS=target,SOURCE_ORIENT=direct)
            if retained:info.update(duplication_metadata(chrom,s,t,a.ref_name,target,direct,
                precise=filler==0 and target0==target1))
            if children:
                ins,deleted=self.edit_lengths(children)
                info.update(DUP_TEMPLATE_LEN=t-s,DUP_SEQUENCE_LEN=qspan[1]-qspan[0],INS_LEN=ins,DEL_LEN=deleted,NET_LEN=q1-q0)
                self.frame_junctions(context,descriptions,a.ref_name,sigma*target0,sigma*target1,a.strand)
                parent=self.make_parent('CPX',a.ref_name,min(sigma*target0,sigma*target1),max(sigma*target0,sigma*target1),qspan,context,children,info)
                if parent:self.insertion_projection(parent,a.ref_name,target,(q0,q1),context,a.strand)
            elif retained:
                self.frame_junctions(context,descriptions,a.ref_name,sigma*target0,sigma*target1,a.strand)
                info.update(NET_LEN=q1-q0,DUP_TEMPLATE_LEN=t-s,DUP_SEQUENCE_LEN=qspan[1]-qspan[0])
                parent=self.make_parent('DUP',chrom,s,t,qspan,context,info=info)
                if parent:self.insertion_projection(parent,a.ref_name,target,(q0,q1),context,a.strand)
            else:
                atom=self.emit_gap(q0,q1,a.ref_name,sigma*target0,sigma*target1,a.strand,context,'backbone')
                if atom:
                    info.update(TARGET_CHR=atom.candidate.chrom,TARGET_POS=atom.candidate.pos)
                    atom.candidate.info.update(info,INS_ORIGIN='TEMPLATED',RESOLUTION='RESOLVED',TEMPLATE_COUNT=1)
                    context.update(resolution='RESOLVED',decomposition='SINGLE')
                    if complete and filler==0 and target0==target1:
                        self.frame_junctions(context,descriptions,a.ref_name,sigma*target0,sigma*target1,a.strand)
                        context['source_join_atoms']=[atom.atom_id]
                else:context['issues'].append('INSERTION_KERNEL_NOT_ACCEPTED')
        elif len(runs)==2 and not children and complete:
            if runs[0]['source'][0]==runs[1]['source'][0] and runs[0]['source'][3]==runs[1]['source'][3]:
                context['issues'].append('SOURCE_GROUPING_AMBIGUOUS');return
            if runs[1]['query'][0]-runs[0]['query'][1]>=SVConfig.MIN_SV_SIZE:
                context['issues'].append('INTERNAL_FILLER_UNEXPLAINED');return
            context['frame']='MULTISOURCE'
            self.frame_junctions(context,descriptions,a.ref_name,sigma*target0,sigma*target1,a.strand)
            signature='|'.join(':'.join(map(str,x['source'])) for x in runs)
            states=['RETAINED' if self.retained(*x['source'][:3],x['query']) else 'UNKNOWN' for x in runs]
            context['source_states']=states
            kinds=['DUP' if state=='RETAINED' else 'INS' for state in states]
            self.make_parent('CPX',a.ref_name,min(sigma*target0,sigma*target1),max(sigma*target0,sigma*target1),(q0,q1),context,
                info=dict(CPX_TYPE='COMPLEX',CPX_COMPONENTS=','.join(kinds),FRAME='MULTISOURCE',INS_ORIGIN='TEMPLATED',
                    TEMPLATE_COUNT=2,SOURCE_STATES=','.join(states),CPX_INTERVALS=','.join(
                        f'{kind}_{x["source"][0]}:{x["source"][1]+1}-{x["source"][2]}' for kind,x in zip(kinds,runs)),
                    INS_LEN=q1-q0,DEL_LEN=0,NET_LEN=q1-q0,SOURCE_SIGNATURE=signature))
        else:context['issues'].append('UNSUPPORTED_SOURCE_PATH')

    def frame_junctions(self,context,segments,chrom,left,right,strand):
        ports=[dict(chrom=chrom,position=left,side='R' if strand=='+' else 'L')]
        pairs=[];last_query=context.get('target_query_span',[segments[0]['query'][0],segments[-1]['query'][1]])[0]
        for segment in segments:
            c,s,t,orientation=segment['source']
            pairs.append((ports[-1],dict(chrom=c,position=s if orientation=='+' else t,side='L' if orientation=='+' else 'R'),[last_query,segment['query'][0]]))
            ports.append(dict(chrom=c,position=t if orientation=='+' else s,side='R' if orientation=='+' else 'L'))
            last_query=segment['query'][1]
        pairs.append((ports[-1],dict(chrom=chrom,position=right,side='L' if strand=='+' else 'R'),[last_query,context.get('target_query_span',[last_query,last_query])[1]]))
        context['junctions']=[junction(self.haplotype,self.name,context['context_id'],a,b,q) for a,b,q in pairs
            if not (a['chrom']==b['chrom'] and a['position']==b['position'] and a['side']!=b['side'])]

    def open_junction(self,context,left,right,members):
        if not members or (left is None and right is None):return
        i,j=(members[-1],right) if left is None else (left,members[0])
        a,b=self.chain[i],self.chain[j]
        if b.qs<a.qe and not self.exact_overlap(a,b,b.qs,a.qe):
            context['issues'].append('QUERY_OWNERSHIP_CONFLICT');return
        p=dict(chrom=a.ref_name,position=exit_port(a),side='R' if a.strand=='+' else 'L')
        q=dict(chrom=b.ref_name,position=entry(b),side='L' if b.strand=='+' else 'R')
        context['junctions']=[junction(self.haplotype,self.name,context['context_id'],p,q,[a.qe,b.qs])]
        edge=context['junctions'][0]
        for suffix,local,remote,mate in [('a',p,q,'b'),('b',q,p,'a')]:
            position=local['position']+(1 if local['side']=='L' else 0)
            pos2=remote['position']+(1 if remote['side']=='L' else 0)
            strands=('+' if local['side']=='R' else '-')+('+' if remote['side']=='L' else '-')
            c=SVCandidate('BND',local['chrom'],max(1,position),0,self.name,
                dict(SOURCE='STRUCTURE',ROLE='JUNCTION',RESOLUTION='PARTIAL',CONTEXT_ID=context['context_id'],
                    CHR2=remote['chrom'],POS2=max(1,pos2),STRANDS=strands,JUNCTION_ID=edge['junction_id']),filter_str='PARTIAL')
            c.id=edge['junction_id']+'.'+suffix;c.mateid=edge['junction_id']+'.'+mate
            self.parents.append(c)

    def annotate_copy(self,atom):
        c=atom.candidate
        if c.type!='INS' or c.info.get('TOPOLOGY') or self.reference is None:return
        mode=SVConfig.COPY_SOURCE_SEARCH
        c.info.setdefault('SOURCE_STATE','UNKNOWN');c.info.setdefault('SOURCE_SEARCH_SCOPE',mode.upper())
        c.info.setdefault('INS_ORIGIN','UNKNOWN')
        if mode=='off':return
        sequence=c.info.get('SEQ','').upper()
        if not sequence or any(x not in 'ACGT' for x in sequence):return
        length=len(sequence);p=c.pos;hits=[]
        if mode in ('adjacent','local'):
            regions=[(p-length,p),(p,p+length)]
            for s,t in regions:
                if s<0:continue
                ref=self.refseq(c.chrom,s,t)
                if ref is None:continue
                if sequence==ref:hits.append((s,t,'DIRECT'))
                if sequence==reverse_complement(ref):hits.append((s,t,'INVERTED'))
            if ((not hits or not any(self.retained(c.chrom,s,t,atom.query) for s,t,_ in hits)) and length<=MAX_LINEAR_BASES and atom.query[1]-atom.query[0]==length and
                    all(a==b for a,b in atom.reference_parts)):
                # A pure insertion may slide through <k identical anchor bases.
                # Verify allele equivalence before testing the normalized port.
                direction=1 if (self.chain[atom.block_index].strand if atom.block_index is not None else self.strand)=='+' else -1
                for distance in range(1,self.k):
                    normalized=[]
                    for delta in (-distance,distance):
                        target=p+delta
                        if target<0:continue
                        pad=self.refseq(c.chrom,min(p,target),max(p,target))
                        if pad is None:continue
                        if delta>0:
                            if sequence[:delta]!=pad:continue
                            shifted=sequence[delta:]+pad
                        else:
                            if sequence[delta:]!=pad:continue
                            shifted=pad+sequence[:delta]
                        span=tuple(q+direction*delta for q in atom.query)
                        if not 0<=span[0]<span[1]<=len(self.sequence):continue
                        observed=self.sequence[span[0]:span[1]].upper()
                        if direction==-1:observed=reverse_complement(observed)
                        if observed!=shifted:continue
                        for s,t in ((target-length,target),(target,target+length)):
                            if s<0:continue
                            ref=self.refseq(c.chrom,s,t)
                            if ref is None:continue
                            if shifted==ref:normalized.append((s,t,'DIRECT',target,span,delta))
                            if shifted==reverse_complement(ref):normalized.append((s,t,'INVERTED',target,span,delta))
                    if normalized:
                        unique=[hit for hit in dict.fromkeys(normalized) if self.retained(c.chrom,hit[0],hit[1],hit[4])]
                        if len(unique)==1:
                            s,t,orientation,target,span,delta=unique[0]
                            atom.source_query=span
                            c.info.update(duplication_metadata(c.chrom,s,t,c.chrom,target,orientation),SOURCE_STATE='RETAINED',
                                INS_ORIGIN='DUPLICATION',TEMPLATE_COUNT=1,DUP_NORMALIZATION='INS_ANCHOR_SHIFT')
                            return
                        if len(unique)>1:
                            c.info.update(SOURCE_STATE='UNKNOWN',SOURCE_ISSUE='SOURCE_AMBIGUOUS');return
        if mode=='local' and not hits:
            s=max(0,p-SVConfig.COPY_SOURCE_RADIUS)
            t=min(self.reference.get_reference_length(c.chrom),p+SVConfig.COPY_SOURCE_RADIUS+length)
            ref=self.refseq(c.chrom,s,t)
            if ref is None:return
            for orientation,text in (('DIRECT',sequence),('INVERTED',reverse_complement(sequence))):
                at=ref.find(text)
                while at>=0 and len(hits)<21:
                    hits.append((s+at,s+at+length,orientation));at=ref.find(text,at+1)
        if len(hits)==1:
            s,t,orientation=hits[0]
            if self.retained(c.chrom,s,t,atom.query):
                c.info.update(duplication_metadata(c.chrom,s,t,c.chrom,p,orientation),
                    SOURCE_STATE='RETAINED',SOURCE_SEARCH_SCOPE=mode.upper(),INS_ORIGIN='DUPLICATION',TEMPLATE_COUNT=1)
        elif len(hits)>1:c.info.update(SOURCE_STATE='UNKNOWN',SOURCE_ISSUE='SOURCE_AMBIGUOUS')

    def promote_duplication(self,atom):
        """Normalize a proven pure CIGAR/gap insertion to a native DUP event.

        The original INS remains unchanged as a reference-coordinate component
        and enters the indel view. Mixed merged representations are annotated
        but are not silently recast as an exact duplication.
        """
        c=atom.candidate
        if atom.parent or c.type!='INS' or c.info.get('TOPOLOGY')!='DUP' or c.info.get('SOURCE_STATE')!='RETAINED':return
        if atom.query[1]-atom.query[0]!=c.sv_len or any(a!=b for a,b in atom.reference_parts):
            c.info.setdefault('SOURCE_ISSUE','MIXED_INSERTION_REPRESENTATION');return
        strand=self.chain[atom.block_index].strand if atom.block_index is not None else self.strand
        source_query=atom.source_query or atom.query
        sequence=self.sequence[source_query[0]:source_query[1]]
        if strand=='-':sequence=reverse_complement(sequence)
        if atom.source_query is None and sequence.upper()!=c.info.get('SEQ','').upper():
            c.info.setdefault('SOURCE_ISSUE','MIXED_INSERTION_REPRESENTATION');return
        index=atom.block_index
        if index is None:
            origin=next((ctx for ctx in self.contexts if ctx['context_id'] in atom.contexts and
                ctx['left_backbone'] is not None and ctx['right_backbone'] is not None),None)
            if origin is None:
                c.info.setdefault('SOURCE_ISSUE','DUP_TARGET_ANCHORS_UNRESOLVED');return
            context=self.context(origin['left_backbone'],origin['right_backbone'],[],'DUP',scope=atom.atom_id)
        else:context=self.context(index,index,[],'DUP',scope=atom.atom_id)
        source=c.info['SOURCE_CHR'];s,t=c.info['SOURCE_START'],c.info['SOURCE_END']
        orientation=strand if c.info['DUP_ORIENT']=='DIRECT' else ('-' if strand=='+' else '+')
        context.update(target_strand=strand,source_intervals=[[source,s,t,orientation]],
            query_span=list(source_query),target_query_span=list(source_query),added_query_span=list(source_query),
            segments=[dict(block_index=index,query=list(source_query),source=[source,s,t,orientation])])
        if atom.source_query is not None:context['insertion_normalization']=dict(original_query=list(atom.query),normalized_query=list(source_query),original_target=c.pos,normalized_target=c.info['TARGET_POS'])
        target=c.info['TARGET_POS']
        self.frame_junctions(context,context['segments'],c.chrom,target,target,strand)
        info={key:value for key,value in c.info.items() if key in ('TOPOLOGY','DUP_MODE','DUP_ORIENT',
            'SUBTYPE','SOURCE_STATE','SOURCE_CHR','SOURCE_START','SOURCE_END','TARGET_CHR','TARGET_POS','SOURCE_SEARCH_SCOPE','DUP_NORMALIZATION')}
        if info['DUP_MODE']=='TANDEM':info['TARGET_POS']=t
        info.update(FRAME='DUP',NET_LEN=c.sv_len,DUP_TEMPLATE_LEN=t-s,DUP_SEQUENCE_LEN=c.sv_len)
        self.make_parent('DUP',source,s,t,source_query,context,info=info,basis=atom,basis_reference=True)

    def run(self):
        if not self.backbone:
            if self.chain:self.branch(None,None,list(range(len(self.chain))))
        else:
            if self.backbone[0]>0:self.branch(None,self.backbone[0],list(range(self.backbone[0])))
            for left,right in zip(self.backbone,self.backbone[1:]):
                if right==left+1:self.ordinary(left,right)
                else:self.branch(left,right,list(range(left+1,right)))
            if self.backbone[-1]<len(self.chain)-1:self.branch(self.backbone[-1],None,list(range(self.backbone[-1]+1,len(self.chain))))
        for i in range(len(self.chain)):self.scan_once(i)
        for atom in self.atoms.values():
            self.annotate_copy(atom)
            if atom.contexts:atom.candidate.info['CONTEXT_ID']=','.join(sorted(atom.contexts))
        for atom in self.atoms.values():self.promote_duplication(atom)
        for atom in self.atoms.values():
            if atom.contexts:atom.candidate.info['CONTEXT_ID']=','.join(sorted(atom.contexts))
        for context in self.contexts:
            context['issues']=sorted(set(context['issues']))
            if context['resolution']!='RESOLVED':
                context['resolution']='PARTIAL' if any(context.get(k) for k in ('atom_ids','segments','junctions')) else 'UNRESOLVED'
        for atom in self.atoms.values():
            if not atom.parent:
                self.events.append(dict(event_id=atom.atom_id,svtype=atom.candidate.type,context_ids=sorted(atom.contexts),
                    query_name=self.name,haplotype=self.haplotype,reference_span=[atom.candidate.chrom,atom.candidate.pos,atom.candidate.end],
                    query_span=list(atom.query),components=[],basis_atom=atom.atom_id,
                    info={k:v for k,v in atom.candidate.info.items() if k!='SEQ'}))
        calls=[atom.candidate for atom in self.atoms.values()]+self.parents
        from asmsv_next.unresolved_regions import annotate_regions
        region_started=time.perf_counter()
        regions=annotate_regions(self,calls+self.projections)
        region_seconds=time.perf_counter()-region_started
        atoms=[dict(atom_id=a.atom_id,type=a.candidate.type,chrom=a.candidate.chrom,pos=a.candidate.pos,
            end=a.candidate.end,length=a.candidate.sv_len,query=list(a.query),reference_consumption=a.reference_parts,
            parts=a.parts,block_index=a.block_index,parent_event=a.parent or None,
            info={k:v for k,v in a.candidate.info.items() if k!='SEQ'},
            contexts=sorted(a.contexts),sequence_sha256=hashlib.sha256(a.candidate.info.get('SEQ','').encode()).hexdigest()) for a in self.atoms.values()]
        return DetectionResult(calls,self.contexts,self.events,atoms,len(self.cache),self.cigar_seconds,self.projections,regions,region_seconds)


def detect_structures(chain,strand,backbone,read_name,sequence,reference=None,haplotype=''):
    return Detector(chain,strand,backbone,read_name,sequence,reference,haplotype).run()
