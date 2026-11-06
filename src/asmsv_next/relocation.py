"""Sequence-compatible source/target links among observed same-haplotype edits.

These relations assert an observed source deletion and matching inserted DNA.
Uniqueness is scoped to observed deletions, never to the entire genome or to a
unique mutational history. Original component calls and atoms are preserved.
"""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from asmsv_next.allele_normalization import shifted_insertion
from asmsv_next.sequence_check import reverse_complement
from asmsv_next.structural_detection import identifier,junction

MAX_LINK_BASES=2000000
MAX_BOUNDARY_SHIFT=15


def rotation_key(sequence):
    """Booth's linear minimal rotation, used only to index possible matches."""
    doubled=sequence+sequence;n=len(sequence);a,b,k=0,1,0
    while a<n and b<n and k<n:
        x,y=doubled[a+k],doubled[b+k]
        if x==y:k+=1;continue
        if x>y:
            a+=k+1
            if a==b:a+=1
        else:
            b+=k+1
            if a==b:b+=1
        k=0
    start=min(a,b)
    return hashlib.sha256(doubled[start:start+n].encode()).hexdigest()


def matching_shifts(sequence,source):
    if sequence==source:return [0]
    if len(sequence)!=len(source):return []
    n=len(sequence);doubled=sequence+sequence;answer=[]
    for shift in range(1,min(MAX_BOUNDARY_SHIFT,n-1)+1):
        if doubled[shift:shift+n]==source:answer.append(shift)
        if doubled[n-shift:2*n-shift]==source:answer.append(-shift)
    return sorted(set(answer),key=lambda x:(abs(x),x))


def deletion_variants(reference,deletion,sequence):
    """Bounded donor shifts proven to leave the full deleted allele equal."""
    yield deletion.pos,deletion.end,sequence,0
    n=len(sequence);chrom=deletion.chrom;s,e=deletion.pos,deletion.end
    for distance in range(1,min(MAX_BOUNDARY_SHIFT,n-1)+1):
        for sign in (-1,1):
            delta=sign*distance
            if s+delta<0:continue
            a,b=(s,s+distance) if sign==1 else (s-distance,s)
            c,d=(e,e+distance) if sign==1 else (e-distance,e)
            first=reference.fetch(chrom,a,b).upper();last=reference.fetch(chrom,c,d).upper()
            if len(first)!=distance or len(last)!=distance or first!=last or any(x not in 'ACGT' for x in first):continue
            rotated=sequence[distance:]+last if sign==1 else first+sequence[:-distance]
            yield s+delta,e+delta,rotated,delta


def link_relocations(calls,haplotype,reference,structure_dir=None):
    index=defaultdict(list);deletions=[]
    for call in calls:
        if call.type!='DEL' or call.filter!='PASS' or call.info.get('ROLE') in ('COMPONENT','PROJECTION'):continue
        if not 0<call.end-call.pos<=MAX_LINK_BASES:continue
        sequence=reference.fetch(call.chrom,call.pos,call.end).upper()
        if len(sequence)!=call.end-call.pos or any(c not in 'ACGT' for c in sequence):continue
        deletions.append((call,sequence))
        index[(len(sequence),rotation_key(sequence))].append((call,sequence,'DIRECT'))
        inverse=reverse_complement(sequence)
        index[(len(sequence),rotation_key(inverse))].append((call,inverse,'INVERTED'))
    relations=[];by_query=defaultdict(list)
    for insertion in calls:
        if insertion.type!='INS' or insertion.filter!='PASS' or insertion.info.get('ROLE') in ('COMPONENT','PROJECTION'):continue
        if insertion.info.get('SOURCE_STATE')=='RETAINED' or insertion.info.get('TOPOLOGY')=='DUP':continue
        sequence=insertion.info.get('SEQ','').upper()
        if not 0<len(sequence)<=MAX_LINK_BASES or any(c not in 'ACGT' for c in sequence):continue
        candidates=[]
        for deletion,original_source,orientation in index.get((len(sequence),rotation_key(sequence)),[]):
            # A same-locus replacement is not evidence of a distant relocation.
            if deletion.chrom==insertion.chrom and deletion.pos<=insertion.pos<=deletion.end:continue
            direct=original_source if orientation=='DIRECT' else reverse_complement(original_source)
            for start,end,dna,donor_shift in deletion_variants(reference,deletion,direct):
                source=dna if orientation=='DIRECT' else reverse_complement(dna)
                for shift in matching_shifts(sequence,source):
                    target=insertion.pos+shift
                    if shifted_insertion(reference,insertion.chrom,insertion.pos,sequence,target)!=source:continue
                    candidates.append((deletion,source,orientation,target,shift,start,end,donor_shift));break
        unique={}
        for match in sorted(candidates,key=lambda m:(abs(m[4])+abs(m[7]),abs(m[4]),m[4],m[7])):
            d,_,o,*_=match
            unique.setdefault((d.chrom,d.pos,d.end,o),match)
        if len(unique)!=1:
            if len(unique)>1:insertion.info.update(SOURCE_ISSUE='OBSERVED_DELETION_SOURCE_AMBIGUOUS',SOURCE_SEARCH_SCOPE='HAPLOTYPE_DELETION_CANDIDATES')
            continue
        deletion,source,orientation,target,shift,start,end,donor_shift=next(iter(unique.values()))
        link=identifier('RELOC',haplotype,insertion.info.get('ATOM_ID',insertion.id),deletion.info.get('ATOM_ID',deletion.id),deletion.chrom,start,end,insertion.chrom,target,orientation)
        insertion.info.update(RELOCATION_ID=link,SOURCE_STATE='REMOVED',INS_ORIGIN='TEMPLATED',
            SOURCE_CHR=deletion.chrom,SOURCE_START=start,SOURCE_END=end,SOURCE_ORIENT=orientation,
            SOURCE_SEARCH_SCOPE='HAPLOTYPE_DELETION_CANDIDATES',TARGET_CHR=insertion.chrom,TARGET_POS=target)
        previous=set(str(deletion.info.get('RELOCATION_ID','')).split(','))-{''}
        deletion.info['RELOCATION_ID']=','.join(sorted(previous|{link}))
        incoming=dict(chrom=deletion.chrom,position=start if orientation=='DIRECT' else end,side='L' if orientation=='DIRECT' else 'R')
        outgoing=dict(chrom=deletion.chrom,position=end if orientation=='DIRECT' else start,side='R' if orientation=='DIRECT' else 'L')
        edges=[junction(haplotype,','.join(sorted(insertion.read_names)),link,dict(chrom=insertion.chrom,position=target,side='R'),incoming),
               junction(haplotype,','.join(sorted(insertion.read_names)),link,outgoing,dict(chrom=insertion.chrom,position=target,side='L'))]
        relation=dict(record_kind='RELOCATION',relation_id=link,haplotype=haplotype,
            scope='HAPLOTYPE_DELETION_CANDIDATES',evidence='EXACT_SEQUENCE_AND_OBSERVED_SOURCE_DELETION',
            historical_mechanism='UNDETERMINED',source=[deletion.chrom,start,end,orientation],
            original_source=[deletion.chrom,deletion.pos,deletion.end],source_boundary_shift=donor_shift,
            target=[insertion.chrom,target],target_boundary_shift=shift,
            insertion_atom=insertion.info.get('ATOM_ID'),deletion_atom=deletion.info.get('ATOM_ID'),
            inserted_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            linked_sequence_sha256=hashlib.sha256(source.encode()).hexdigest(),junctions=edges)
        relations.append(relation)
        for query in insertion.read_names:by_query[query].append(relation)
    if structure_dir:
        for path in Path(structure_dir).glob(haplotype+'_*.json'):
            bundle=json.loads(path.read_text());links=by_query.get(bundle['query'],[])
            if not links:continue
            bundle['relocations']=links
            # Derived links supersede a flattened insertion adjacency only in
            # the adjacency projection; the original atoms remain in JSON.
            path.write_text(json.dumps(bundle,ensure_ascii=False,allow_nan=False)+'\n')
    return relations
