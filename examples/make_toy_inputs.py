#!/usr/bin/env python3
"""Write a small installation fixture with explicit, synthetic alignments.

Hap1 contains a 200 bp deletion and a 60 bp insertion; Hap2 is reference.
This is an installation smoke example, not an accuracy benchmark or a mapper.
"""
import argparse
from pathlib import Path
import random

import pysam


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    rng=random.Random(20261101)
    reference=''.join(rng.choice('ACGT') for _ in range(10000))
    inserted=''.join(rng.choice('ACGT') for _ in range(60))
    fasta=args.output/'reference.fa'
    fasta.write_text('>chrToy\n'+reference+'\n')
    pysam.faidx(str(fasta))
    header=dict(HD={'VN':'1.6','SO':'queryname'},SQ=[{'SN':'chrToy','LN':len(reference)}],
                PG=[{'ID':'minimap2','PN':'minimap2','VN':'2.28','CL':'minimap2 -ax asm5 --eqx reference.fa assembly.fa'}])
    for hap in (1,2):
        sequence=reference[:1000]+reference[1200:2200]+inserted+reference[2200:3200] if hap==1 else reference[:3200]
        name=f'hap{hap}_toy'
        assembly=args.output/f'hap{hap}.fa'
        assembly.write_text('>'+name+'\n'+sequence+'\n');pysam.faidx(str(assembly))
        with pysam.AlignmentFile(str(args.output/f'hap{hap}.bam'),'wb',header=header) as bam:
            row=pysam.AlignedSegment(bam.header)
            row.query_name=name;row.query_sequence=sequence;row.flag=0
            row.reference_id=0;row.reference_start=0;row.mapping_quality=60
            row.cigartuples=[(7,1000),(2,200),(7,1000),(1,60),(7,1000)] if hap==1 else [(7,3200)]
            row.set_tag('NM',260 if hap==1 else 0);bam.write(row)
    print(f'Created toy inputs in {args.output.resolve()}')


if __name__=='__main__':main()
