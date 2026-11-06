"""Whole-contig BAM groups and optional lazy original-assembly sequence access."""
from dataclasses import dataclass, replace, field
from functools import lru_cache
from itertools import groupby
from pathlib import Path
from multiprocessing import shared_memory
import pysam
from asmsv_next.alignment_tags import source_from_header
from asmsv_next.config import SVConfig
from asmsv_next.presets import canonical_aligner, mode_from_command


@lru_cache(maxsize=4)
def fasta_handle(path):return pysam.FastaFile(path)


@dataclass(frozen=True)
class QuerySequence:
    name: str
    length: int
    sequence: str = None
    fasta_path: str = None
    shared_name: str = None
    _shared_handle: object = field(default=None, init=False, repr=False, compare=False)

    def __len__(self):return self.length

    def __getitem__(self,item):
        if isinstance(item,int):
            if item<0:item+=self.length
            if not 0<=item<self.length:raise IndexError(item)
            start,end=item,item+1
        else:
            start,end,step=item.indices(self.length)
            if step!=1:raise ValueError('Query provider requires unit-stride slices')
        if end<=start:return ''
        if self.sequence is not None:return self.sequence[start:end]
        if self.shared_name is not None:
            handle=self._shared_handle
            if handle is None:
                handle=shared_memory.SharedMemory(name=self.shared_name)
                object.__setattr__(self,'_shared_handle',handle)
            return bytes(handle.buf[start:end]).decode('ascii')
        return fasta_handle(self.fasta_path).fetch(self.name,start,end)

    def close_shared(self):
        if self._shared_handle is not None:
            self._shared_handle.close()
            object.__setattr__(self,'_shared_handle',None)


def inspected_source(header,aligner,mode=None,program_id=None):
    source=source_from_header(header,program_id)
    if source.aligner!='unknown':
        observed=canonical_aligner(source.aligner)
        if observed!=aligner:raise ValueError(f'Declared aligner {aligner} conflicts with BAM source {source.aligner}')
        observed_mode=mode_from_command(source.command_line)
        if mode and observed_mode and mode!=observed_mode:
            raise ValueError(f'Declared source preset {mode} conflicts with BAM preset {observed_mode}')
        source=replace(source,mode=observed_mode)
    return source


def inspect_headers(paths,reference,aligner):
    ref={name:reference.get_reference_length(name) for name in reference.references}
    sources=[]
    for path in paths:
        with pysam.AlignmentFile(str(path),'rb') as bam:
            observed={x['SN']:x['LN'] for x in bam.header.to_dict().get('SQ',[])}
            if not observed:raise ValueError('Missing BAM reference dictionary: '+str(path))
            if any(ref.get(name)!=size for name,size in observed.items()):
                raise ValueError('BAM and FASTA reference names/lengths disagree: '+str(path))
            sources.append(inspected_source(bam.header.to_dict(),aligner))
    return sources


def _query_provider(reads,name,length,query_fasta):
    if query_fasta:
        fasta=fasta_handle(query_fasta)
        if name not in fasta.references or fasta.get_reference_length(name)!=length:
            raise ValueError(f'Original assembly FASTA missing {name} or contig length differs')
        samples=[r for r in reads if r.query_length and r.query_sequence is not None]
        if samples:
            r=min(samples,key=lambda x:x.query_length)
            sequence=r.get_forward_sequence()
            clips=r.cigartuples or []
            left_h=clips[0][1] if clips and clips[0][0]==5 else 0
            right_h=clips[-1][1] if clips and clips[-1][0]==5 else 0
            start=right_h if r.is_reverse else left_h
            end=start+len(sequence)
            width=min(64,len(sequence))
            if (sequence[:width].upper()!=fasta.fetch(name,start,start+width).upper() or
                    sequence[-width:].upper()!=fasta.fetch(name,end-width,end).upper()):
                raise ValueError('Assembly FASTA/BAM sequence endpoint sample differs: '+name)
        return QuerySequence(name,length,fasta_path=query_fasta)
    complete=[r for r in reads if r.query_length==length and r.query_sequence is not None]
    if not complete:return None
    r=min(complete,key=lambda x:(bool(x.flag&0x900),-x.mapping_quality,x.flag,x.reference_id,x.reference_start))
    return QuerySequence(name,length,sequence=r.get_forward_sequence())


def bam_groups(path,haplotype,query_fasta=None,sequence_transform=None):
    """Input must be queryname sorted; every group retains its original coordinates."""
    p=SVConfig.CALLER_PRESET
    query_fasta=str(Path(query_fasta).resolve()) if query_fasta else None
    source_id=str(Path(path).resolve())
    with pysam.AlignmentFile(str(path),'rb') as bam:
        header=bam.header.to_dict();cache={};seen=set()
        for name,group in groupby(bam,key=lambda r:r.query_name):
            if name in seen:raise ValueError('BAM query groups are not contiguous; name-sort the BAM: '+name)
            seen.add(name)
            reads=list(group)
            if any(r.has_tag('CG') and isinstance(r.get_tag('CG'),str) for r in reads):
                raise ValueError(f'{name}: text CG:Z is unsupported; retain real CIGAR, avoid VACmap --L')
            mapped=[r for r in reads if not r.is_unmapped and r.cigartuples]
            if not mapped:continue
            if not any(SVConfig.MIN_MAPQ<=r.mapping_quality<255 for r in mapped):continue
            lengths={sum(n for op,n in r.cigartuples if op in (0,1,4,5,7,8)) for r in mapped}
            if len(lengths)!=1:raise ValueError('BAM full query coordinates disagree: '+name)
            length=lengths.pop();sequence=_query_provider(reads,name,length,query_fasta)
            if sequence is None:
                raise ValueError(f'{name}: no complete SEQ; provide --query-fasta1/2 for hard-clipped assembly output')
            records=[]
            for r in mapped:
                pg=r.get_tag('PG') if r.has_tag('PG') else None
                if pg not in cache:cache[pg]=inspected_source(header,p.aligner,p.mode,pg)
                tags={k:r.get_tag(k) for k in ('AS','NM','MD','cs','nn','PG') if r.has_tag(k)}
                records.append(dict(q_name=name,ref_name=r.reference_name,ref_start=r.reference_start,
                    ref_end=r.reference_end,cigar=r.cigartuples,mapq=r.mapping_quality,
                    is_reverse=r.is_reverse,flag=r.flag,query_length=length,raw_as=tags.get('AS'),
                    tags=tags,source=cache[pg],source_id=source_id,haplotype=haplotype))
            if sequence_transform is not None:sequence=sequence_transform(sequence)
            # No native BAM record is needed after its immutable payload has
            # been built. In particular, release cached whole-contig SEQ before
            # the generator suspends while the worker runs.
            del reads,mapped,r
            yield name,records,sequence,length,haplotype
