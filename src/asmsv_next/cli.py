"""AsmSV: assembly-based haplotype-resolved structural variant calling."""
import argparse
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import uuid
import importlib.util
import shutil
import tempfile
import threading
import time
import pysam
from asmsv_next import __version__
from asmsv_next.config import SVConfig
from asmsv_next.constants import fixed_policy
from asmsv_next.genotype import cluster_and_genotype,cluster_structural_calls
from asmsv_next.input_io import inspect_headers
from asmsv_next.pipeline import run_haplotype
from asmsv_next.presets import canonical_aligner,mode_from_command,resolve_preset
from asmsv_next.runtime import arm_parent_death,raise_termination
from asmsv_next.structures import BedMaskFilter
from asmsv_next.vcf import write_vcf

logger=logging.getLogger(__name__)


def arguments(argv=None):
    p=argparse.ArgumentParser(prog='asmSV',description=__doc__)
    p.add_argument('--version',action='version',version='AsmSV '+__version__)
    p.add_argument('--aligner',required=True,type=canonical_aligner,choices=('minimap2','winnowmap','vacmap'))
    p.add_argument('--aligner-preset',help='Upstream mode for provenance/validation; inferred from @PG or asm5/asm defaults; does not select caller scores')
    p.add_argument('--preset-file',help='Optional schema-5 JSON override of the fixed common AQ/DP template; upstream scoring flags are never inherited')
    p.add_argument('--bam1',required=True);p.add_argument('--bam2',required=True);p.add_argument('--ref',required=True)
    p.add_argument('--query-fasta1');p.add_argument('--query-fasta2')
    p.add_argument('--sorted-bam1',help='Existing BAM grouped by queryname; never modified')
    p.add_argument('--sorted-bam2',help='Existing BAM grouped by queryname; never modified')
    p.add_argument('--output','-o',default='results.vcf');p.add_argument('--threads','-t',type=int,default=4)
    p.add_argument('--work-dir',help='Override the queryname BAM cache directory; by default each cache is stored beside its input BAM')
    p.add_argument('--sequence-transport',choices=('auto','copy'),default='auto',help='Share large BAM-derived query sequences when /dev/shm has space; copy uses the original process transport')
    p.add_argument('--scratch-dir',help='Temporary BAM sorting and structural evidence directory; defaults to the system temporary directory. Final JSONL evidence is retained with the VCF')
    p.add_argument('--trace-dir',help='Optional per-contig AQ, DP and candidate diagnostics without sequences')
    p.add_argument('--min-mapq',type=int,default=10);p.add_argument('--min-sv',type=int,default=30)
    p.add_argument('--min-inv','--min-inv-size',dest='min_inv',type=int,help='Minimum reported standalone INV size; defaults to --min-sv; does not change CPX frame interpretation')
    p.add_argument('--max-indel-size',type=int,default=-1,help='Inclusive reporting size limit for INS/DEL records; larger records retain their type with OVERSIZE; -1 disables')
    p.add_argument('--max-ref-gap',type=int,default=150000,help='Reference-gap limit for heuristic fallback only; exact pure indels and CIGAR evidence are separate; not a TRA classification rule')
    p.add_argument('--max-sv',type=int,default=1500000,help='Retained legacy duplication OVERSIZE threshold, not a universal size filter')
    p.add_argument('--mask',help='Exclude candidates overlapping this BED; not an annotation-only input')
    p.add_argument('--debug',action='store_true')
    p.add_argument('--replacement-realign',choices=('off','fallback'),default='off',help='Optional bounded WFA after unresolved k-mer recovery; requires pywfa')
    p.add_argument('--seed-max-bases',type=int,default=20000)
    p.add_argument('--seed-max-blocks',type=int,default=1000)
    p.add_argument('--dup-source-search','--copy-source-search',dest='copy_source_search',choices=('off','adjacent','local'),default='adjacent')
    p.add_argument('--dup-source-radius','--copy-source-radius',dest='copy_source_radius',type=int,default=100000,help='Optional local source-start search radius; does not define DUP')
    p.add_argument('--junction-vcf',action='store_true',help='Write an additional read-only BND projection; does not add main events')
    p.add_argument('--tandem-target-indels',action='store_true',help='Also write a reference-indel view with exactly equivalent tandem INS at the declared copy target; benchmark representation must be specified')
    p.add_argument('--vcf-info',choices=('minimal','full'),default='minimal',help='Biological annotations by default; full also exports diagnostic fields')
    p.add_argument('--sample', default='SAMPLE', help='VCF sample name, e.g. HG002')
    p.add_argument('--genotyping', choices=('assembly','evidence','presence'), default='assembly', help='Structured detector: assembly retains dosage estimates with FT confidence flags; evidence makes unsupported counterparts missing; presence reproduces the previous policy')
    a=p.parse_args(argv)
    if not a.sample or any(c.isspace() for c in a.sample):p.error('Sample name must be nonempty and contain no whitespace')
    if a.threads<1 or not 0<=a.min_mapq<255 or a.min_sv<1:
        p.error('Invalid worker count, minimum size or MAPQ')
    if a.max_sv!=-1 and a.max_sv<a.min_sv:p.error('max-sv is smaller than min-sv')
    if a.min_inv is not None and a.min_inv<1:p.error('min-inv must be positive')
    if a.max_indel_size!=-1 and a.max_indel_size<a.min_sv:p.error('max-indel-size is smaller than min-sv')
    if a.max_ref_gap<1:p.error('max-ref-gap must be positive')
    if a.seed_max_bases<1 or a.seed_max_blocks<1 or a.copy_source_radius<0:p.error('Invalid sequence-check budget')
    if a.replacement_realign=='fallback' and importlib.util.find_spec('pywfa') is None:p.error('--replacement-realign fallback requires the optional pywfa package')
    return a


def save_json(path,data):
    path=Path(path);tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    try:
        tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
    finally:
        if tmp.exists():tmp.unlink()


def identity(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)


def file_digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


@contextmanager
def run_stage(manifest,record,name):
    """Expose long native-library calls without relying on their stderr output."""
    manifest['stage']=name;save_json(record,manifest)
    started=time.monotonic();stopped=threading.Event()
    def heartbeat():
        while not stopped.wait(30):
            logger.info('%s: still running (%.1f s)',name,time.monotonic()-started)
    thread=threading.Thread(target=heartbeat,daemon=True)
    logger.info('Starting %s',name);thread.start()
    success=False
    try:
        yield
        success=True
    finally:
        stopped.set();thread.join()
        elapsed=time.monotonic()-started
        manifest.setdefault('stage_seconds',{})[name]=elapsed
        save_json(record,manifest)
        logger.info('%s %s in %.2f s','Finished' if success else 'Stopped',name,elapsed)


def sorted_bam(path,provided,work,threads,scratch_dir=None):
    if provided:
        if not Path(provided).is_file():raise FileNotFoundError(provided)
        pysam.quickcheck(str(provided))
        logger.info('Using supplied queryname BAM: %s',provided)
        return str(Path(provided).resolve())
    proof=identity(path);key=hashlib.sha256(json.dumps(proof,sort_keys=True).encode()).hexdigest()[:20]
    source=Path(proof['path'])
    if work is None:
        out=source.with_name(source.stem+'.asmsv.qname.bam')
        receipt=Path(str(out)+'.json')
    else:
        work=Path(work).expanduser().resolve()
        work.mkdir(parents=True,exist_ok=True)
        # Preserve the pre-1.3.1 explicit --work-dir layout and valid receipts.
        out=work/(key+'.qname.bam');receipt=work/(key+'.json')
    if out.resolve()==source:raise ValueError('Sorting cache must not overwrite the input BAM')
    if out.exists() and receipt.exists():
        try:
            saved=json.loads(receipt.read_text())
            if saved==proof or (saved.get('schema')=='asmsv-sort-cache-v1' and
                    saved.get('source')==proof and saved.get('cache')==identity(out)):
                pysam.quickcheck(str(out))
                with pysam.AlignmentFile(str(out),'rb') as stream:
                    if stream.header.to_dict().get('HD',{}).get('SO')!='queryname':
                        raise ValueError('Cached BAM does not declare queryname sorting')
                logger.info('Reusing queryname BAM cache: %s',out);return str(out)
            logger.info('Sorting cache is stale; rebuilding %s',out)
        except (OSError,ValueError,AttributeError,pysam.SamtoolsError) as error:
            logger.warning('Sorting cache validation failed for %s: %s; rebuilding',out,error)
    elif out.exists():
        logger.info('Sorting cache has no source receipt; rebuilding %s',out)
    tmp=out.with_name(out.name+'.'+uuid.uuid4().hex+'.tmp.bam')
    try:
        try:tmp.touch(exist_ok=False)
        except PermissionError as error:
            raise PermissionError(f'Cannot create sorting cache beside {path}; use --work-dir with a writable directory or --sorted-bam1/2') from error
        # Spill files can be several times larger than the input. Keep them off
        # an output-local network filesystem, while publishing the cache atomically.
        with tempfile.TemporaryDirectory(prefix='asmsv-sort-',dir=scratch_dir) as sorting:
            logger.info('Name-sorting %s; temporary files: %s; cache: %s',path,sorting,out)
            # Native pysam sorting can hold the GIL, suppressing progress and
            # delaying KeyboardInterrupt. A child keeps the controller responsive;
            # subprocess.run kills/reaps it on interruption before scratch cleanup.
            subprocess.run([sys.executable,'-c',
                'import sys,pysam; pysam.sort(*sys.argv[1:],catch_stdout=False)',
                '-n','-@',str(threads),'-m','256M','-T',str(Path(sorting)/'chunk'),
                '-o',str(tmp),str(path)],check=True)
        pysam.quickcheck(str(tmp))
        if identity(path)!=proof:raise ValueError('Input BAM changed during name sorting; retry with a stable input: '+str(path))
        os.replace(tmp,out)
        save_json(receipt,dict(schema='asmsv-sort-cache-v1',source=proof,cache=identity(out)))
    finally:
        if tmp.exists():tmp.unlink()
    return str(out)


def run(argv=None):
    started=time.monotonic()
    a=arguments(argv);out=Path(a.output).expanduser().resolve()
    if out.exists() and out.is_dir():out=out/'results.vcf'
    out.parent.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.DEBUG if a.debug else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(str(out)+'.log')],force=True)
    if os.environ.get('ASMSV_OWNER_PID'):arm_parent_death(int(os.environ['ASMSV_OWNER_PID']))
    signal.signal(signal.SIGTERM,raise_termination);signal.signal(signal.SIGHUP,raise_termination)
    manifest=dict(status='RUNNING',software_version=__version__,started_at=datetime.now(timezone.utc).isoformat(),argv=list(argv or sys.argv[1:]),
        implementation=str(Path(__file__).resolve().parent),inputs={},csv_policy='finite-structure-v2',
        annotation_schema='asmsv-sv-annotation-v3',vcf_info=a.vcf_info,
        runtime=dict(python=sys.version,pysam=pysam.__version__,samtools=pysam.__samtools_version__,
            libdeflate=bool(getattr(pysam.config,'HAVE_LIBDEFLATE',False))),stage_seconds={},
        genotyping=a.genotyping,sample=a.sample)
    manifest['implementation_sha256']={path.name:file_digest(path) for path in sorted(Path(__file__).resolve().parent.glob('*.py'))}
    record=Path(str(out)+'.run.json')
    scratch=None;sequences=None
    try:
        for key in ('bam1','bam2','ref','query_fasta1','query_fasta2','mask'):
            value=getattr(a,key)
            if value:manifest['inputs'][key]=identity(value)
        save_json(record,manifest)
        logger.info('AsmSV %s; pysam %s; samtools %s; libdeflate=%s',__version__,
            pysam.__version__,pysam.__samtools_version__,manifest['runtime']['libdeflate'])
        if not manifest['runtime']['libdeflate']:
            logger.warning('This pysam build lacks libdeflate; BAM sorting/compression can be substantially slower. See docs/USAGE.md')
        with pysam.FastaFile(a.ref) as ref:
            sources=inspect_headers((a.bam1,a.bam2),ref,a.aligner)
            modes={mode_from_command(s.command_line) for s in sources if mode_from_command(s.command_line)}
            if len(modes)>1:raise ValueError('Hap1/Hap2 upstream presets differ; use matching source modes for one run')
            mode=a.aligner_preset or (next(iter(modes)) if modes else ('asm' if a.aligner=='vacmap' else 'asm5'))
            if modes and modes!={mode}:raise ValueError('User upstream preset conflicts with BAM @PG')
            preset=resolve_preset(a.aligner,mode,preset_file=a.preset_file)
            SVConfig.CALLER_PRESET=preset;SVConfig.MIN_MAPQ=a.min_mapq
            SVConfig.MIN_SV_SIZE=a.min_sv;SVConfig.MAX_SV_SIZE=a.max_sv
            SVConfig.MIN_INV_SIZE=a.min_inv;SVConfig.MAX_INDEL_SIZE=a.max_indel_size
            SVConfig.MIN_TRA_SIZE=a.max_ref_gap
            SVConfig.DETECTOR='structured';SVConfig.REFERENCE_FASTA=str(Path(a.ref).resolve())
            SVConfig.SEED_MAX_BASES=a.seed_max_bases;SVConfig.SEED_MAX_BLOCKS=a.seed_max_blocks
            SVConfig.REPLACEMENT_REALIGN=a.replacement_realign
            SVConfig.VCF_INFO=a.vcf_info
            SVConfig.GENOTYPING=a.genotyping
            SVConfig.SAMPLE_NAME=a.sample
            SVConfig.SEQUENCE_TRANSPORT=a.sequence_transport
            SVConfig.COPY_SOURCE_SEARCH=a.copy_source_search;SVConfig.COPY_SOURCE_RADIUS=a.copy_source_radius
            manifest['detector_config']={key:value for key,value in vars(SVConfig).items()
                if key.isupper() and key!='CALLER_PRESET'}
            manifest.update(source_headers=[asdict(s) for s in sources],preset=preset.to_dict(),
                derived_parameters=preset.derived_parameters(),fixed_policy=fixed_policy(),
                preset_sha256=preset.fingerprint,scoring_sha256=preset.scoring_fingerprint,
                scoring_policy='explicit-json' if a.preset_file else 'fixed-common',
                upstream_scoring_inherited=False,
                mode_origin='user' if a.aligner_preset else 'header' if modes else 'declared-tool-default')
            save_json(str(out)+'.preset.json',preset.to_dict());save_json(record,manifest)
            logger.info('Resolved %s; score unit=%s; preset SHA256=%s',preset.id,preset.score_unit,preset.fingerprint)
            work=Path(a.work_dir).expanduser().resolve() if a.work_dir else None
            if a.scratch_dir:Path(a.scratch_dir).mkdir(parents=True,exist_ok=True)
            scratch=tempfile.TemporaryDirectory(prefix='asmsv-structures-',dir=a.scratch_dir)
            structure_dir=Path(scratch.name);SVConfig.STRUCTURE_DIR=str(structure_dir)
            manifest['structure_dir']=str(structure_dir)
            manifest['structure_scratch_retained']=False
            manifest['detector_config']['STRUCTURE_DIR']=str(structure_dir)
            paths=[]
            for label,path,provided in (('hap1',a.bam1,a.sorted_bam1),('hap2',a.bam2,a.sorted_bam2)):
                with run_stage(manifest,record,'prepare_'+label):
                    paths.append(sorted_bam(path,provided,work,a.threads,a.scratch_dir))
                    manifest['sorted_bams']=list(paths)
            inspect_headers(paths,ref,a.aligner)
            trace=str(Path(a.trace_dir).resolve()) if a.trace_dir else None
            if trace:Path(trace).mkdir(parents=True,exist_ok=True)
            mask=BedMaskFilter(a.mask)
            with run_stage(manifest,record,'call_hap1'):
                h1,s1=run_haplotype(paths[0],'hap1',a.threads,mask,a.query_fasta1,trace)
            with run_stage(manifest,record,'call_hap2'):
                h2,s2=run_haplotype(paths[1],'hap2',a.threads,mask,a.query_fasta2,trace)
            from asmsv_next.relocation import link_relocations
            relations=link_relocations(h1,'hap1',ref,SVConfig.STRUCTURE_DIR)+link_relocations(h2,'hap2',ref,SVConfig.STRUCTURE_DIR)
            manifest['relocation_summary']=dict(relations=len(relations),scope='HAPLOTYPE_DELETION_CANDIDATES')
            writing_started=time.monotonic()
            manifest['stage']='genotype_and_write';save_json(record,manifest)
            logger.info('Starting cross-haplotype genotyping and output')
            from asmsv_next.allele_normalization import canonical_insertions
            (variants, aliases) = cluster_structural_calls(canonical_insertions(h1,ref), canonical_insertions(h2,ref))
            for c in variants:
                if c.info.get('PARENT_EVENT'):
                    c.info['PARENT_EVENT']=','.join(sorted({aliases.get(k,k) for k in c.info['PARENT_EVENT'].split(',')}))
            contigs={name:ref.get_reference_length(name) for name in ref.references}
            primary=[c for c in variants if c.info.get('ROLE')!='COMPONENT']
            from asmsv_next.genotype_evidence import GenotypeEvidence,summarize_genotypes
            evidence=None
            if a.genotyping!='presence':
                bundles=(dict(haplotype=s['haplotype'],genotyping=s.pop('_genotyping')) for s in s1+s2)
                evidence=GenotypeEvidence(bundles,flank=preset.aq.min_anchor)
            genotype_summary={}
            record_maps=[]
            region_links=defaultdict(list)
            from asmsv_next.sequence_output import SequenceStore
            sequences=SequenceStore(str(out)+'.sequences.fa')
            def output_view(calls,path,view):
                from asmsv_next.structures import apply_indel_size_policy
                apply_indel_size_policy(calls)
                if evidence is not None:evidence.apply(calls,policy=a.genotyping)
                for candidate in calls:sequences.annotate(candidate)
                genotype_summary[view]=summarize_genotypes(calls)
                mapping=Path(SVConfig.STRUCTURE_DIR,view+'.records.jsonl')
                with mapping.open('w',buffering=4*1024*1024) as records:
                    def save_record(data):
                        records.write(json.dumps(data,ensure_ascii=False,allow_nan=False)+'\n')
                        for region_id in str(data['info'].get('REGION_ID','')).split(','):
                            if region_id:region_links[region_id].append((view,data['record_id']))
                    write_vcf(calls,str(path),contigs,ref,view=view,record_sink=save_record)
                record_maps.append(mapping)
            output_view(primary,out,'main')
            components=Path(str(out)+'.components.vcf');indels=Path(str(out)+'.indels.vcf');graph=Path(str(out)+'.events.jsonl')
            output_view([c for c in variants if c.info.get('ROLE')=='COMPONENT'],components,'components')
            # A source-coordinate edit in an added copy is not a reference-locus DEL.
            def reference_indels(calls):
                return [c for c in calls if c.type in ('INS','DEL') and not (
                    c.info.get('ROLE')=='COMPONENT' and c.info.get('ALLELE_CONTEXT') in ('ADDED_DUPLICATE','DUPLICATION_JUNCTION'))]
            from asmsv_next.allele_normalization import canonical_tandem_indels
            parents={c.id:c for c in h1+h2 if c.type=='DUP'}
            normalized_haps=[canonical_insertions(reference_indels(calls),ref,indel_view=True) for calls in (h1,h2)]
            projected=cluster_and_genotype(
                *normalized_haps)
            for c in projected:
                if c.info.get('PARENT_EVENT'):
                    c.info['PARENT_EVENT']=','.join(sorted({aliases.get(k,k) for k in c.info['PARENT_EVENT'].split(',')}))
            output_view(projected,indels,'indels')
            target_indels=None
            if a.tandem_target_indels:
                target_indels=Path(str(out)+'.target_indels.vcf')
                target_calls=cluster_and_genotype(*(canonical_tandem_indels(calls,parents,ref) for calls in normalized_haps))
                for candidate in target_calls:
                    if candidate.info.get('PARENT_EVENT'):
                        candidate.info['PARENT_EVENT']=','.join(sorted({aliases.get(k,k) for k in candidate.info['PARENT_EVENT'].split(',')}))
                output_view(target_calls,target_indels,'target_indels')
            region_spool=Path(SVConfig.STRUCTURE_DIR,'regions.spool.jsonl')
            with graph.open('w',buffering=4*1024*1024) as stream,region_spool.open('w',buffering=4*1024*1024) as region_stream:
                junction_haps={'hap1':[],'hap2':[]}
                stream.write(json.dumps(dict(record_kind='EVENT_ALIASES',aliases=aliases),ensure_ascii=False)+'\n')
                for relation in relations:stream.write(json.dumps(relation,ensure_ascii=False,allow_nan=False)+'\n')
                for mapping in record_maps:
                    with mapping.open() as records:shutil.copyfileobj(records,stream)
                for path in sorted(Path(SVConfig.STRUCTURE_DIR).glob('*.json')):
                    data=json.loads(path.read_text());data['record_kind']='CONTIG'
                    stream.write(json.dumps(data,ensure_ascii=False,allow_nan=False)+'\n')
                    for region in data.get('regions',[]):
                        region_stream.write(json.dumps(region,ensure_ascii=False,allow_nan=False)+'\n')
                    if a.junction_vcf:
                        from asmsv_next.output_views import junction_candidates
                        junction_haps[data['haplotype']].extend(junction_candidates(data,aliases))
            manifest['views']=dict(components=str(components),indels=str(indels),events=str(graph))
            if target_indels is not None:manifest['views']['target_indels']=str(target_indels)
            if a.junction_vcf:
                junction_path=Path(str(out)+'.junctions.vcf')
                output_view(cluster_and_genotype(junction_haps['hap1'],junction_haps['hap2']),junction_path,'junctions')
                with graph.open('a',buffering=4*1024*1024) as stream,record_maps[-1].open() as records:shutil.copyfileobj(records,stream)
                manifest['views']['junctions']=str(junction_path)
            from asmsv_next.region_output import write_regions
            region_report=write_regions(SVConfig.STRUCTURE_DIR,record_maps,out,aliases,contigs,mask_applied=bool(a.mask),record_links=region_links,region_source=region_spool)
            manifest['views'].update(regions=region_report['jsonl'],regions_bed=region_report['bed'])
            manifest['region_summary']=region_report
            sequences.finish()
            manifest['views']['sequences']=str(sequences.path)
            if sequences.seen:manifest['views']['sequence_index']=str(sequences.path)+'.fai'
            manifest['sequence_summary']=dict(unique_payloads=len(sequences.seen),bases=sum(sequences.seen.values()))
            genotype_report=Path(str(out)+'.genotypes.json')
            save_json(genotype_report,dict(schema='asmsv-genotypes-v1',sample=a.sample,policy=a.genotyping,views=genotype_summary))
            manifest['genotype_summary']=genotype_summary
            manifest['genotype_report']=str(genotype_report)
            manifest['genotype_report_sha256']=file_digest(genotype_report)
            manifest['view_sha256']={k:file_digest(v) for k,v in manifest['views'].items()}
            with pysam.VariantFile(str(out)) as stream:rows=sum(1 for _ in stream)
            manifest.update(status='COMPLETE',contigs=s1+s2,output_rows=rows,
                output_sha256=file_digest(out),sorted_bams=paths)
            manifest['stage_seconds']['genotype_and_write']=time.monotonic()-writing_started
            manifest['stage']='complete'
            logger.info('Completed: %s; %d VCF records',out,rows)
    except BaseException as error:
        manifest.update(status='FAILED',error=f'{type(error).__name__}: {error}')
        logger.exception('Caller failed')
        raise
    finally:
        manifest['elapsed_seconds']=time.monotonic()-started
        manifest['finished_at']=datetime.now(timezone.utc).isoformat();save_json(record,manifest)
        if sequences is not None:sequences.abort()
        if scratch is not None:scratch.cleanup()
    return 0


def main():
    try:return run()
    except (Exception,KeyboardInterrupt):return 1


if __name__=='__main__':sys.exit(main())
