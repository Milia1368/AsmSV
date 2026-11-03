"""Bounded whole-contig tasks: AQ, global DP, backbone and structural detection."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import asdict
import hashlib
import json
import multiprocessing
import os
import time
from pathlib import Path
from asmsv_next.alignment import AlignedSegment
from asmsv_next.config import SVConfig
from asmsv_next.dp import solve_path,select_backbone,edge_terms
from asmsv_next.input_io import bam_groups,fasta_handle
from asmsv_next.runtime import abort_pool,worker_lifetime
from asmsv_next.structures import filter_candidates_by_mask
from asmsv_next.sequence_transport import SequenceTransport


def initialize_worker(config,parent_pid):
    worker_lifetime(parent_pid)
    for key,value in config.items():setattr(SVConfig,key,value)


def describe(block):
    return dict(query=[block.qs,block.qe],reference=[block.ref_name,block.ref_start,block.ref_end,block.strand],
        native_AS=block.raw_as,mapq=block.mapq,mapq_factor=block.mapq_factor,
        AQ_sequence=block.aq_sequence,AQ_effective=block.aq,E_effective=block.as_score,
        summary=asdict(block.summary),rejection_reasons=block.rejection_reasons,
        quality_flags=block.quality_flags,reward_rescued=block.reward_rescued,
        alignment_id=block.score_root.raw.alignment_id,
        source_query=[block.score_root.raw.qs,block.score_root.raw.qe],
        source=asdict(block.score_root.raw.source))


def call_contig(name,records,sequence,length,haplotype,trace_dir=None):
    try:
        return _call_contig_body(name,records,sequence,length,haplotype,trace_dir)
    finally:
        sequence.close_shared()


def _call_contig_body(name,records,sequence,length,haplotype,trace_dir=None):
    t0=time.perf_counter()
    blocks=[AlignedSegment(r) for r in records];valid=[b for b in blocks if b.is_valid]
    input_count,eligible_count=len(blocks),len(valid)
    reasons=Counter(reason for b in blocks for reason in b.rejection_reasons)
    t1=time.perf_counter()
    result=solve_path(valid,length);chain=list(result.blocks)
    t2=time.perf_counter()
    indices,strand=select_backbone(chain,length)
    t3=time.perf_counter()
    if not trace_dir:
        # AQ/DP have consumed these nodes. The selected chain owns the roots
        # needed downstream; unselected roots need not live through SV parsing.
        del blocks,valid
    structural=None
    from asmsv_next.structural_detection import detect_structures
    reference=fasta_handle(SVConfig.REFERENCE_FASTA) if SVConfig.REFERENCE_FASTA else None
    structural=detect_structures(chain,strand,indices,name,sequence,reference,haplotype)
    calls=structural.calls+structural.projections;cigar_seconds=structural.cigar_seconds
    t4=time.perf_counter()-cigar_seconds
    stats=dict(haplotype=haplotype,query=name,length=length,input_blocks=input_count,eligible_blocks=eligible_count,
        chain_blocks=len(chain),backbone_blocks=len(indices),branch_blocks=len(chain)-len(indices),
        candidates=len(calls),dp_score=result.score,evaluated_predecessors=result.evaluated_predecessors,
        tied=result.tied,rejections=dict(reasons),query_sequence='FASTA' if sequence.fasta_path else 'BAM')
    stats['sequence_backend']='SHARED' if sequence.shared_name else 'FASTA' if sequence.fasta_path else 'COPY'
    stats['stage_seconds']=dict(aq=t1-t0,global_dp=t2-t1,backbone_dp=t3-t2,
        split_detection=t4-t3-structural.region_seconds,cigar_detection=cigar_seconds,
        region_annotation=structural.region_seconds)
    if structural:
        stats.update(cigar_scans=structural.cigar_scans,structure_contexts=len(structural.contexts),
            resolved_cpx=sum(c.type=='CPX' for c in calls),
            unresolved_contexts=sum(c['resolution']!='RESOLVED' for c in structural.contexts),
            unresolved_regions=len(structural.regions))
        if SVConfig.STRUCTURE_DIR:
            key=hashlib.sha256((haplotype+'\0'+name).encode()).hexdigest()[:20]
            from asmsv_next.genotype_evidence import contig_evidence
            data=dict(schema='asmsv-structure-v2',query=name,haplotype=haplotype,
                blocks=[dict(index=i,alignment_id=b.score_root.raw.alignment_id,query=[b.qs,b.qe],
                    reference=[b.ref_name,b.ref_start,b.ref_end,b.strand]) for i,b in enumerate(chain)],
                backbone_indices=sorted(indices),contexts=structural.contexts,events=structural.events,atoms=structural.atoms,
                region_schema='asmsv-regions-v1',regions=structural.regions,
                genotyping=contig_evidence(chain,indices,name))
            if SVConfig.GENOTYPING!='presence':stats['_genotyping']=data['genotyping']
            target=Path(SVConfig.STRUCTURE_DIR,f'{haplotype}_{key}.json')
            tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(data,ensure_ascii=False,allow_nan=False)+'\n');os.replace(tmp,target)
    if trace_dir:
        key=hashlib.sha256((haplotype+'\0'+name).encode()).hexdigest()[:20]
        data=dict(stats={k:v for k,v in stats.items() if k!='_genotyping'},preset=SVConfig.CALLER_PRESET.to_dict(),blocks=[describe(b) for b in blocks],
            chain=[describe(b) for b in chain],backbone_indices=sorted(indices),
            transitions=[dict(left=a.score_root.raw.alignment_id,right=b.score_root.raw.alignment_id,
                terms=edge_terms(a,b)) for a,b in zip(chain,chain[1:])],
            candidates=[dict(type=c.type,chrom=c.chrom,pos=c.pos,end=c.end,length=c.sv_len,
                query_interval=[c.q_start,c.q_end] if c.q_start or c.q_end else None,
                reference_gap=c.reference_gap,
                source=c.info.get('SOURCE'),subtype=c.info.get('SUBTYPE')) for c in calls])
        Path(trace_dir,f'{haplotype}_{key}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    return calls,stats


def run_haplotype(path,label,workers,mask=None,query_fasta=None,trace_dir=None):
    config={k:v for k,v in vars(SVConfig).items() if k.isupper()}
    all_calls=[];stats=[];pending={}
    transport=SequenceTransport(SVConfig.SEQUENCE_TRANSPORT)
    pool=ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
        initializer=initialize_worker,initargs=(config,os.getpid()))
    def consume(done):
        for future in done:
            shared_name=pending.pop(future)
            try:
                calls,summary=future.result();all_calls.extend(calls);stats.append(summary)
            finally:transport.release(shared_name)
    try:
        for payload in bam_groups(path,label,query_fasta,sequence_transform=transport.prepare):
            future=pool.submit(call_contig,*payload,trace_dir)
            pending[future]=payload[2].shared_name
            del payload
            if len(pending)>=workers:
                done,_=wait(pending,return_when=FIRST_COMPLETED);consume(done)
        while pending:
            done,_=wait(pending,return_when=FIRST_COMPLETED);consume(done)
    except BaseException:
        abort_pool(pool,pending)
        raise
    finally:
        pool.shutdown(wait=True)
        transport.close()
    all_calls.sort(key=lambda c:(tuple(sorted(c.read_names)),c.chrom,c.pos,c.type,c.sv_len,c.end,c.id))
    if mask and mask.is_active:all_calls=filter_candidates_by_mask(all_calls,mask)
    return all_calls,sorted(stats,key=lambda x:x['query'])
