"""AsmSV: assembly-based haplotype-resolved structural variant calling."""
import argparse
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import sys
import uuid
import importlib.util
import shutil
import tempfile
import pysam
from asmsv_next import __version__
from asmsv_next.config import SVConfig
from asmsv_next.constants import fixed_policy
from asmsv_next.genotype import cluster_and_genotype,cluster_structural_calls
from asmsv_next.input_io import inspect_headers
from asmsv_next.pipeline import run_haplotype
from asmsv_next.presets import canonical_aligner,mode_from_command,resolve_preset,source_costs
from asmsv_next.runtime import arm_parent_death,raise_termination
from asmsv_next.structures import BedMaskFilter
from asmsv_next.vcf import write_vcf

logger=logging.getLogger(__name__)


def arguments(argv=None):
    p=argparse.ArgumentParser(prog='asmSV',description=__doc__)
    p.add_argument('--version',action='version',version='AsmSV '+__version__)
    p.add_argument('--aligner',required=True,type=canonical_aligner,choices=('minimap2','winnowmap','vacmap'))
    p.add_argument('--aligner-preset',help='Upstream mode; inferred from unambiguous @PG or asm5/asm defaults')
    p.add_argument('--preset-file',help='Optional schema-5 JSON override of the remaining AQ/DP parameters')
    p.add_argument('--bam1',required=True);p.add_argument('--bam2',required=True);p.add_argument('--ref',required=True)
    p.add_argument('--query-fasta1');p.add_argument('--query-fasta2')
    p.add_argument('--sorted-bam1',help='Existing BAM grouped by queryname; never modified')
    p.add_argument('--sorted-bam2',help='Existing BAM grouped by queryname; never modified')
    p.add_argument('--output','-o',default='results.vcf');p.add_argument('--threads','-t',type=int,default=4)
    p.add_argument('--work-dir',help='Own queryname sort cache; defaults to output directory/cache')
    p.add_argument('--sequence-transport',choices=('auto','copy'),default='auto',help='Share large BAM-derived query sequences when /dev/shm has space; copy uses the original process transport')
    p.add_argument('--scratch-dir',help='Temporary structural evidence directory; defaults to the system temporary directory. Final JSONL evidence is retained with the VCF')
    p.add_argument('--trace-dir',help='Optional per-contig AQ, DP and candidate diagnostics without sequences')
    p.add_argument('--min-mapq',type=int,default=10);p.add_argument('--min-sv',type=int,default=30)
    p.add_argument('--max-sv',type=int,default=1500000,help='Retained legacy duplication OVERSIZE threshold, not a universal size filter')
    p.add_argument('--mask',help='Exclude candidates overlapping this BED; not an annotation-only input')
    p.add_argument('--debug',action='store_true')
    p.add_argument('--replacement-realign',choices=('off','fallback'),default='off',help='Optional bounded WFA after unresolved k-mer recovery; requires pywfa')
    p.add_argument('--seed-max-bases',type=int,default=20000)
    p.add_argument('--seed-max-blocks',type=int,default=1000)
    p.add_argument('--dup-source-search','--copy-source-search',dest='copy_source_search',choices=('off','adjacent','local'),default='adjacent')
    p.add_argument('--dup-source-radius','--copy-source-radius',dest='copy_source_radius',type=int,default=100000,help='Optional local source-start search radius; does not define DUP')
    p.add_argument('--junction-vcf',action='store_true',help='Write an additional read-only BND projection; does not add main events')
    p.add_argument('--vcf-info',choices=('minimal','full'),default='minimal',help='Biological annotations by default; full also exports diagnostic fields')
    p.add_argument('--sample', default='SAMPLE', help='VCF sample name, e.g. HG002')
    p.add_argument('--genotyping', choices=('assembly','evidence','presence'), default='assembly', help='Structured detector: assembly retains dosage estimates with FT confidence flags; evidence makes unsupported counterparts missing; presence reproduces the previous policy')
    a=p.parse_args(argv)
    if not a.sample or any(c.isspace() for c in a.sample):p.error('Sample name must be nonempty and contain no whitespace')
    if a.threads<1 or not 0<=a.min_mapq<255 or a.min_sv<1:
        p.error('Invalid worker count, minimum size or MAPQ')
    if a.max_sv!=-1 and a.max_sv<a.min_sv:p.error('max-sv is smaller than min-sv')
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


def sorted_bam(path,provided,work,threads):
    if provided:
        if not Path(provided).is_file():raise FileNotFoundError(provided)
        return str(Path(provided).resolve())
    proof=identity(path);key=hashlib.sha256(json.dumps(proof,sort_keys=True).encode()).hexdigest()[:20]
    out=work/(key+'.qname.bam');receipt=work/(key+'.json')
    if out.exists() and receipt.exists() and json.loads(receipt.read_text())==proof:
        pysam.quickcheck(str(out));return str(out)
    tmp=work/(key+'.'+uuid.uuid4().hex+'.bam')
    try:
        pysam.sort('-n','-@',str(threads),'-m','256M','-T',str(tmp)+'.sort','-o',str(tmp),str(path))
        pysam.quickcheck(str(tmp));os.replace(tmp,out);save_json(receipt,proof)
    finally:
        if tmp.exists():tmp.unlink()
    return str(out)


def run(argv=None):
    a=arguments(argv);out=Path(a.output).expanduser().resolve()
    if out.exists() and out.is_dir():out=out/'results.vcf'
    out.parent.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.DEBUG if a.debug else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(str(out)+'.log')],force=True)
    if os.environ.get('ASMSV_OWNER_PID'):arm_parent_death(int(os.environ['ASMSV_OWNER_PID']))
    signal.signal(signal.SIGTERM,raise_termination);signal.signal(signal.SIGHUP,raise_termination)
    manifest=dict(status='RUNNING',started_at=datetime.now(timezone.utc).isoformat(),argv=list(argv or sys.argv[1:]),
        implementation=str(Path(__file__).resolve().parent),inputs={},csv_policy='finite-structure-v2',
        annotation_schema='asmsv-sv-annotation-v3',vcf_info=a.vcf_info,
        genotyping=a.genotyping,sample=a.sample)
    manifest['implementation_sha256']={path.name:file_digest(path) for path in sorted(Path(__file__).resolve().parent.glob('*.py'))}
    record=Path(str(out)+'.run.json')
    scratch=None
    try:
        for key in ('bam1','bam2','ref','query_fasta1','query_fasta2','mask'):
            value=getattr(a,key)
            if value:manifest['inputs'][key]=identity(value)
        save_json(record,manifest)
        with pysam.FastaFile(a.ref) as ref:
            sources=inspect_headers((a.bam1,a.bam2),ref,a.aligner)
            modes={mode_from_command(s.command_line) for s in sources if mode_from_command(s.command_line)}
            if len(modes)>1:raise ValueError('Hap1/Hap2 upstream presets differ; use matching source modes for one run')
            mode=a.aligner_preset or (next(iter(modes)) if modes else ('asm' if a.aligner=='vacmap' else 'asm5'))
            if modes and modes!={mode}:raise ValueError('User upstream preset conflicts with BAM @PG')
            overrides=[s.overrides for s in sources]
            costs=[source_costs(a.aligner,mode,v) for v in overrides]
            if costs[0]!=costs[1] and not a.preset_file:
                raise ValueError('Hap1/Hap2 source scoring differs; provide an explicit common --preset-file')
            preset=resolve_preset(a.aligner,mode,overrides[0],a.preset_file)
            SVConfig.CALLER_PRESET=preset;SVConfig.MIN_MAPQ=a.min_mapq
            SVConfig.MIN_SV_SIZE=a.min_sv;SVConfig.MAX_SV_SIZE=a.max_sv
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
                preset_sha256=preset.fingerprint,mode_origin='user' if a.aligner_preset else 'header' if modes else 'declared-tool-default')
            save_json(str(out)+'.preset.json',preset.to_dict());save_json(record,manifest)
            logger.info('Resolved %s; score unit=%s; preset SHA256=%s',preset.id,preset.score_unit,preset.fingerprint)
            work=Path(a.work_dir).resolve() if a.work_dir else out.parent/'cache'
            work.mkdir(parents=True,exist_ok=True)
            if a.scratch_dir:Path(a.scratch_dir).mkdir(parents=True,exist_ok=True)
            scratch=tempfile.TemporaryDirectory(prefix='asmsv-structures-',dir=a.scratch_dir)
            structure_dir=Path(scratch.name);SVConfig.STRUCTURE_DIR=str(structure_dir)
            manifest['structure_dir']=str(structure_dir)
            manifest['structure_scratch_retained']=False
            manifest['detector_config']['STRUCTURE_DIR']=str(structure_dir)
            paths=[sorted_bam(a.bam1,a.sorted_bam1,work,a.threads),sorted_bam(a.bam2,a.sorted_bam2,work,a.threads)]
            inspect_headers(paths,ref,a.aligner)
            trace=str(Path(a.trace_dir).resolve()) if a.trace_dir else None
            if trace:Path(trace).mkdir(parents=True,exist_ok=True)
            mask=BedMaskFilter(a.mask)
            h1,s1=run_haplotype(paths[0],'hap1',a.threads,mask,a.query_fasta1,trace)
            h2,s2=run_haplotype(paths[1],'hap2',a.threads,mask,a.query_fasta2,trace)
            (variants, aliases) = cluster_structural_calls(h1, h2)
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
            def output_view(calls,path,view):
                if evidence is not None:evidence.apply(calls,policy=a.genotyping)
                genotype_summary[view]=summarize_genotypes(calls)
                mapping=Path(SVConfig.STRUCTURE_DIR,view+'.records.jsonl')
                with mapping.open('w',buffering=4*1024*1024) as records:
                    def save_record(data):records.write(json.dumps(data,ensure_ascii=False,allow_nan=False)+'\n')
                    write_vcf(calls,str(path),contigs,ref,view=view,record_sink=save_record)
                record_maps.append(mapping)
            output_view(primary,out,'main')
            components=Path(str(out)+'.components.vcf');indels=Path(str(out)+'.indels.vcf');graph=Path(str(out)+'.events.jsonl')
            output_view([c for c in variants if c.info.get('ROLE')=='COMPONENT'],components,'components')
            # A source-coordinate edit in an added copy is not a reference-locus DEL.
            def reference_indels(calls):
                return [c for c in calls if c.type in ('INS','DEL') and not (
                    c.info.get('ROLE')=='COMPONENT' and c.info.get('ALLELE_CONTEXT') in ('ADDED_DUPLICATE','DUPLICATION_JUNCTION'))]
            projected=cluster_and_genotype(reference_indels(h1),reference_indels(h2))
            for c in projected:
                if c.info.get('PARENT_EVENT'):
                    c.info['PARENT_EVENT']=','.join(sorted({aliases.get(k,k) for k in c.info['PARENT_EVENT'].split(',')}))
            output_view(projected,indels,'indels')
            with graph.open('w',buffering=4*1024*1024) as stream:
                junction_haps={'hap1':[],'hap2':[]}
                stream.write(json.dumps(dict(record_kind='EVENT_ALIASES',aliases=aliases),ensure_ascii=False)+'\n')
                for mapping in record_maps:
                    with mapping.open() as records:shutil.copyfileobj(records,stream)
                for path in sorted(Path(SVConfig.STRUCTURE_DIR).glob('*.json')):
                    data=json.loads(path.read_text());data['record_kind']='CONTIG'
                    stream.write(json.dumps(data,ensure_ascii=False,allow_nan=False)+'\n')
                    if a.junction_vcf:
                        from asmsv_next.output_views import junction_candidates
                        junction_haps[data['haplotype']].extend(junction_candidates(data,aliases))
            manifest['views']=dict(components=str(components),indels=str(indels),events=str(graph))
            if a.junction_vcf:
                junction_path=Path(str(out)+'.junctions.vcf')
                output_view(cluster_and_genotype(junction_haps['hap1'],junction_haps['hap2']),junction_path,'junctions')
                with graph.open('a',buffering=4*1024*1024) as stream,record_maps[-1].open() as records:shutil.copyfileobj(records,stream)
                manifest['views']['junctions']=str(junction_path)
            genotype_report=Path(str(out)+'.genotypes.json')
            save_json(genotype_report,dict(schema='asmsv-genotypes-v1',sample=a.sample,policy=a.genotyping,views=genotype_summary))
            manifest['genotype_summary']=genotype_summary
            manifest['genotype_report']=str(genotype_report)
            manifest['genotype_report_sha256']=file_digest(genotype_report)
            manifest['view_sha256']={k:file_digest(v) for k,v in manifest['views'].items()}
            with pysam.VariantFile(str(out)) as stream:rows=sum(1 for _ in stream)
            manifest.update(status='COMPLETE',contigs=s1+s2,output_rows=rows,
                output_sha256=file_digest(out),sorted_bams=paths)
            logger.info('Completed: %s; %d VCF records',out,rows)
    except BaseException as error:
        manifest.update(status='FAILED',error=f'{type(error).__name__}: {error}')
        logger.exception('Caller failed')
        raise
    finally:
        manifest['finished_at']=datetime.now(timezone.utc).isoformat();save_json(record,manifest)
        if scratch is not None:scratch.cleanup()
    return 0


def main():
    try:return run()
    except (Exception,KeyboardInterrupt):return 1


if __name__=='__main__':sys.exit(main())
