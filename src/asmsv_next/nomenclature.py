"""Observable structure names; GATK-SV names are used only for exact templates.

CPX/DELINS are declared caller extensions, not additional core VCF ALT types.
COPY is never an emitted variant class. Unknown nomenclature is independent of
whether the selected structure is fully resolved.
"""


def duplication_metadata(source_chrom,start,end,target_chrom,target,orientation,
                         precise=True):
    adjacent=source_chrom==target_chrom and target in (start,end)
    if adjacent and orientation=='DIRECT' and precise:
        mode,subtype='TANDEM','DUP:TANDEM'
    elif adjacent:
        mode='ADJACENT'
        subtype=('piDUP_FR' if target==end else 'piDUP_RF') if orientation=='INVERTED' and precise else 'DUP'
    else:mode,subtype='DISPERSED','dDUP'
    return dict(TOPOLOGY='DUP',DUP_MODE=mode,DUP_ORIENT=orientation,SUBTYPE=subtype,
        SOURCE_CHR=source_chrom,SOURCE_START=start,SOURCE_END=end,
        TARGET_CHR=target_chrom,TARGET_POS=target)


def complex_metadata(frame,sources,children):
    """Do not equate internal edits with named flanking-deletion templates."""
    components=['INV' if frame=='INV' else 'DUP'] if frame in ('INV','DUP') else []
    intervals=[]
    for source in sources:
        chrom,start,end=source[:3]
        if components and end>start:intervals.append(f'{components[0]}_{chrom}:{start+1}-{end}')
    for atom in sorted(children,key=lambda a:(a.candidate.chrom,a.candidate.pos,a.candidate.end,a.candidate.type)):
        kind=atom.candidate.type
        components.extend(['DEL','INS'] if kind=='DELINS' else [kind])
        if kind in ('DEL','DELINS'):
            for start,end in atom.reference_parts:
                if end>start:intervals.append(f'DEL_{atom.candidate.chrom}:{start+1}-{end}')
    subtype='COMPLEX'
    if frame=='INV' and len(sources)==1 and children and all(a.candidate.type=='DEL' for a in children):
        chrom,start,end=sources[0][:3];left=right=False;internal=False
        for atom in children:
            if atom.candidate.chrom!=chrom:internal=True;continue
            for a,b in atom.reference_parts:
                if b<=start:left=True
                elif a>=end:right=True
                else:internal=True
        if not internal:
            subtype='delINVdel' if left and right else 'delINV' if left else 'INVdel' if right else 'COMPLEX'
    info=dict(CPX_TYPE=subtype,CPX_COMPONENTS=','.join(components))
    if intervals:info['CPX_INTERVALS']=','.join(intervals)
    return info
