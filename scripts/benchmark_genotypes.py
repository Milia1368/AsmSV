#!/usr/bin/env python3
"""One-to-one INS/DEL genotype benchmark, independent of ALT-count reuse.

The historical detection benchmark uses --pick ac. Genotypes are additionally
evaluated using --pick single, so dosage does not control match multiplicity.
FT-filtered estimates stay in the benchmark; their retained fraction is reported.
"""
import argparse
import json
from pathlib import Path

import pysam

from evaluate_benchmark import BCF,BENCH_PY,REF,TRUTH,BED,command,digest,save
from evaluate_genotypes import evaluate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('vcf',type=Path,help='Reference-change INS/DEL VCF')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ref',required=True);p.add_argument('--truth',required=True);p.add_argument('--bed',required=True)
    p.add_argument('--resume',action='store_true',help='Reuse completed results only after verifying input and result hashes')
    a=p.parse_args();out=a.output.resolve()
    args=['--pick','single','--passonly','-d','-r','2000','-C','5000','-s','50','-S','30',
          '--sizemax','50000','-p','0.7','-P','0.7','-O','0.0']
    signature=dict(input_vcf=str(a.vcf.resolve()),vcf_sha256=digest(a.vcf),truth_sha256=digest(a.truth),
                   bed_sha256=digest(a.bed),reference_sha256=digest(a.ref),args=args)
    if out.exists() and any(out.iterdir()):
        if not a.resume:p.error('Use a fresh output directory or --resume')
        if not (out/'evaluation.json').is_file() or not (out/'inputs.json').is_file():p.error('Incomplete benchmark; use a fresh output directory')
        previous=json.loads((out/'evaluation.json').read_text())
        if json.loads((out/'inputs.json').read_text())!=signature:p.error('Benchmark inputs changed')
        if any(digest(out/'bench'/name)!=expected for name,expected in previous['inputs'].items()):p.error('Benchmark results changed')
        print('Reused benchmark after verifying input and result hashes: '+str(out));return
    out.mkdir(parents=True,exist_ok=True)
    plain=out/'calls.vcf';compressed=out/'calls.vcf.gz';bench=out/'bench'
    with pysam.VariantFile(str(a.vcf)) as stream,pysam.VariantFile(str(plain),'w',header=stream.header) as target:
        for row in stream:
            if row.info.get('SVTYPE') in ('INS','DEL') and set(row.filter)=={'PASS'}:target.write(row)
    command(out,'sort',[BCF,'sort',plain,'-Oz','-o',compressed]);command(out,'index',[BCF,'index','-f','-t',compressed])
    save(out/'inputs.json',signature)
    command(out,'bench',[BENCH_PY,'-m','truvari','bench','-b',a.truth,'-c',compressed,'-o',bench,
        *args,'-f',a.ref,'--includebed',a.bed])
    result,base,comp=evaluate(bench)
    with (out/'pairs.jsonl').open('w') as stream:
        for side,rows in [('base',base),('comp',comp)]:
            for row in rows:stream.write(json.dumps(dict(side=side,**row))+'\n')
    save(out/'evaluation.json',result)
    print(json.dumps({k:result[k] for k in ('detection','genotype_aware','base','comp')},indent=2))


if __name__=='__main__':main()
