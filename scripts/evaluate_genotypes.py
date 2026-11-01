#!/usr/bin/env python3
"""Genotype concordance from completed Truvari matches, including missing alleles.

GTMatch in Truvari is an alternate-allele-count difference, so it is insufficient
to distinguish 0/1 from a partial 1/. call. Explicit GT tuples are compared here.
No haplotype-order or phase accuracy is claimed by this evaluator.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import pysam


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def gt_name(gt):
    return '/'.join('.' if a is None else str(a) for a in sorted(gt,key=lambda a:(a is None, a if a is not None else -1)))


def complete(gt):
    return bool(gt) and all(a is not None for a in gt)


def rows(path):
    result=[]
    with pysam.VariantFile(str(path)) as stream:
        for row in stream:
            gt=tuple(row.samples[0].get('GT') or ())
            length=row.info.get('SVLEN',row.stop-row.pos)
            if isinstance(length,tuple):length=length[0]
            mid=row.info.get('MatchId',())
            if isinstance(mid,str):mid=(mid,)
            result.append(dict(chrom=row.chrom,pos=row.pos,id=row.id,kind=row.info.get('SVTYPE','UNKNOWN'),
                size=abs(length),gt=gt,match=tuple(mid),phased=row.samples[0].phased,
                ft=row.samples[0].get('FT') if 'FT' in row.format else None,
                ps=row.samples[0].get('PS') if 'PS' in row.format else None))
    return result


def comparison(base,comp):
    full=complete(base['gt']) and complete(comp['gt'])
    same=full and gt_name(base['gt'])==gt_name(comp['gt'])
    return 'concordant' if same else 'discordant' if full else 'partial_or_missing'


def paired_summary(own,other,base_side):
    own_key,other_key=(0,1) if base_side else (1,0)
    lookup={}
    for row in other:
        if len(row['match'])!=2:raise ValueError('Expected Truvari MatchId=(base-key,comp-key)')
        key=row['match'][other_key]
        if key in lookup:raise ValueError('Truvari side has duplicate variant key: '+key)
        lookup[key]=row
    counts=Counter();matrix=defaultdict(Counter);strata=defaultdict(Counter)
    qualified=Counter();filters=Counter();sizes=defaultdict(Counter)
    pairs=[]
    for row in own:
        if len(row['match'])!=2:raise ValueError('Expected two-part MatchId')
        key=row['match'][other_key]
        if key not in lookup:raise ValueError('Missing matched Truvari record: '+key)
        partner=lookup[key]
        base,comp=(row,partner) if base_side else (partner,row)
        label=comparison(base,comp)
        counts[label]+=1;matrix[gt_name(base['gt'])][gt_name(comp['gt'])]+=1
        strata[base['kind']][label]+=1
        if comp['ft'] in (None,'.','PASS'):qualified[label]+=1
        else:filters[str(comp['ft'])]+=1
        size=base['size']
        bucket=next((name for bound,name in ((100,'<100'),(300,'100-299'),(1000,'300-999'),(10000,'1k-9999')) if size<bound),'>=10k')
        sizes[bucket][label]+=1
        pairs.append(dict(chrom=base['chrom'],truth_pos=base['pos'],call_pos=comp['pos'],call_id=comp['id'],
            type=base['kind'],truth_gt=gt_name(base['gt']),call_gt=gt_name(comp['gt']),genotype_filter=comp['ft'],result=label,match=list(row['match'])))
    count=len(own);full=counts['concordant']+counts['discordant']
    qtotal=sum(qualified.values())
    return dict(records=count,**{k:counts[k] for k in ('concordant','discordant','partial_or_missing')},
        full_genotype_fraction=full/count if count else None,
        concordance_fully_called=counts['concordant']/full if full else None,
        concordance_all_matched=counts['concordant']/count if count else None,
        unfiltered_genotypes=dict(records=qtotal,**{k:qualified[k] for k in ('concordant','discordant','partial_or_missing')},
            fraction=qtotal/count if count else None,concordance=qualified['concordant']/qtotal if qtotal else None),
        genotype_filters=dict(filters),by_size={k:dict(v) for k,v in sorted(sizes.items())},
        matrix={k:dict(v) for k,v in sorted(matrix.items())},by_type={k:dict(v) for k,v in sorted(strata.items())}),pairs


def evaluate(bench):
    summary=json.loads((bench/'summary.json').read_text())
    base=rows(bench/'tp-base.vcf.gz');comp=rows(bench/'tp-comp.vcf.gz')
    base_result,base_pairs=paired_summary(base,comp,True)
    comp_result,comp_pairs=paired_summary(comp,base,False)
    p=comp_result['concordant']/(summary['TP-comp']+summary['FP']) if summary['TP-comp']+summary['FP'] else 0
    r=base_result['concordant']/(summary['TP-base']+summary['FN']) if summary['TP-base']+summary['FN'] else 0
    matrices={}
    for label in ('fp','fn'):
        entries=rows(bench/(label+'.vcf.gz'))
        matrices[label]=dict(Counter(gt_name(row['gt']) for row in entries))
    result=dict(schema='asmsv-genotype-evaluation-v1',benchmark=str(bench.resolve()),
        match_policy=json.loads((bench/'params.json').read_text()) if (bench/'params.json').exists() else {},
        detection={k:summary[k] for k in ('TP-base','TP-comp','FP','FN','precision','recall','f1')},
        genotype_aware=dict(precision=p,recall=r,f1=2*p*r/(p+r) if p+r else 0),
        base=base_result,comp=comp_result,unmatched_genotypes=matrices,
        truvari_ac_concordance=summary.get('gt_concordance'),
        caveats=['Genotypes are compared without allele order; this is not a phasing assessment.',
                 'Assembly-policy genotypes with FT flags are dosage estimates; unfiltered concordance and the retained fraction are also reported.',
                 'Partial genotypes are failures in all-matched and genotype-aware metrics, and excluded only from the explicitly conditional fully-called concordance.',
                 'Truvari --pick ac can use different base/comp multiplicities; both denominators are reported.',
                 'Truth is INS/DEL; native DUP/INV/CPX genotype accuracy requires additional truth.'],
        inputs={label+'.vcf.gz':digest(bench/(label+'.vcf.gz')) for label in ('tp-base','tp-comp','fp','fn')})
    return result,base_pairs,comp_pairs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('benchmark',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--compare',type=Path,help='Previous genotype evaluation JSON')
    a=p.parse_args();result,base_pairs,comp_pairs=evaluate(a.benchmark)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    detail=a.output.with_suffix('.pairs.jsonl')
    with detail.open('w') as stream:
        for label,pairs in [('base',base_pairs),('comp',comp_pairs)]:
            for row in pairs:stream.write(json.dumps(dict(side=label,**row),ensure_ascii=False)+'\n')
    result['pairs_file']=str(detail.resolve())
    if a.compare:
        previous=json.loads(a.compare.read_text());result['comparison']={}
        for key in ('precision','recall','f1'):
            result['comparison'][key+'_delta_pp']=100*(result['genotype_aware'][key]-previous['genotype_aware'][key])
        for key in ('concordance_fully_called','concordance_all_matched','full_genotype_fraction'):
            old=previous['base'][key];new=result['base'][key]
            result['comparison'][key+'_delta_pp']=100*(new-old) if old is not None and new is not None else None
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('detection','genotype_aware','base','comp','comparison') if k in result},indent=2))


if __name__=='__main__':main()
