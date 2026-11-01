#!/usr/bin/env python3
"""Benchmark a completed caller output; retain commands, hashes and truth identities."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import sys
import time

import pysam
from evaluate_genotypes import evaluate as evaluate_genotypes

BENCH_PY = sys.executable
BCF = shutil.which('bcftools') or 'bcftools'
REF = TRUTH = BED = None
BENCH_ARGS = ['--pick', 'ac', '--passonly', '-d', '-r', '2000', '-C', '5000']


def save(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def command(directory, name, argv):
    argv = list(map(str, argv))
    before = time.monotonic()
    with (directory / (name + '.stdout.log')).open('w') as out, (directory / (name + '.stderr.log')).open('w') as err:
        result = subprocess.run(argv, stdout=out, stderr=err, check=False)
    save(directory / (name + '.command.json'), dict(argv=argv, returncode=result.returncode,
                                                  seconds=time.monotonic() - before))
    if result.returncode:
        raise RuntimeError(f'{name} failed: {directory / (name + ".stderr.log")}')


def variant_key(row):
    return json.dumps([row.chrom, row.pos, row.id, row.ref, row.alts], separators=(',', ':'))


def summarize(bench):
    raw = json.loads((bench / 'summary.json').read_text())
    result = {k: raw[k] for k in ('TP-base', 'TP-comp', 'FP', 'FN', 'precision', 'recall', 'f1')}
    genotypes,_,_ = evaluate_genotypes(bench)
    result['genotype_aware'] = genotypes['genotype_aware']
    result['genotype_concordance'] = {side: {key: genotypes[side][key] for key in (
        'records','concordant','discordant','partial_or_missing','full_genotype_fraction',
        'concordance_fully_called','concordance_all_matched','matrix','by_type',
        'unfiltered_genotypes','genotype_filters','by_size')} for side in ('base','comp')}
    result['truvari_ac_concordance'] = raw.get('gt_concordance')
    result['strata'] = {}
    for label in ('tp-base', 'tp-comp', 'fp', 'fn'):
        types, sizes, chromosomes = Counter(), Counter(), Counter()
        with pysam.VariantFile(str(bench / (label + '.vcf.gz'))) as stream:
            for row in stream:
                kind = row.info.get('SVTYPE', 'UNKNOWN')
                size = row.info.get('SVLEN', row.stop - row.pos)
                size = abs(size[0] if isinstance(size, tuple) else size)
                bucket = next((name for bound, name in ((50, '<50'), (100, '50-99'), (300, '100-299'),
                             (1000, '300-999'), (10000, '1k-9999')) if size < bound), '>=10k')
                types[kind] += 1; sizes[bucket] += 1; chromosomes[row.chrom] += 1
        result['strata'][label] = dict(types=dict(types), sizes=dict(sizes), chromosomes=dict(chromosomes))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf', type=Path)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--truth', required=True)
    parser.add_argument('--bed', required=True)
    parser.add_argument('--compare', type=Path, help='Previously benchmarked caller directory')
    args = parser.parse_args()
    directory = args.vcf.resolve().parent
    manifest_path = Path(str(args.vcf) + '.run.json')
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest['status'] != 'COMPLETE' or manifest['output_sha256'] != digest(args.vcf):
            raise ValueError('Caller manifest is incomplete or VCF has changed')
    signature = dict(vcf_sha256=digest(args.vcf), truth_sha256=digest(args.truth),projection_policy='exclude_csv_cpx_and_delins_v2',
                     bed_sha256=digest(args.bed), reference_sha256=digest(args.ref), args=BENCH_ARGS)
    def view(name):
        target = directory / (name + '.vcf.gz')
        bench = directory / ('bench_' + name)
        proof = directory / (name + '.input.json')
        if bench.exists():
            if not proof.exists() or json.loads(proof.read_text()) != signature:
                raise ValueError('Existing benchmark has different inputs: ' + str(bench))
            if (bench / 'summary.json').exists():
                return name, summarize(bench)
            raise ValueError('Incomplete benchmark: ' + str(bench))
        plain = directory / (name + '.vcf')
        with pysam.VariantFile(str(args.vcf)) as source, pysam.VariantFile(str(plain), 'w', header=source.header) as out:
            for row in source:
                kind = row.info.get('SVTYPE')
                if set(row.filter) == {'PASS'} and kind not in ('CSV','CPX','DELINS') and (name == 'primary' or kind in ('INS', 'DEL')):
                    out.write(row)
        command(directory, name + '_sort', [BCF, 'sort', plain, '-Oz', '-o', target])
        command(directory, name + '_index', [BCF, 'index', '-f', '-t', target])
        save(proof, signature)
        command(directory, name + '_bench', [BENCH_PY, '-m', 'truvari', 'bench', '-b', args.truth,
                '-c', target, '-o', bench, *BENCH_ARGS, '-f', args.ref, '--includebed', args.bed])
        return name, summarize(bench)
    with ThreadPoolExecutor(max_workers=2) as pool:
        metrics = dict(pool.map(view, ('primary', 'common_insdel')))
    save(directory / 'metrics.json', metrics)
    if args.compare:
        changes = {}
        for name in metrics:
            sets = []
            for folder in (args.compare, directory):
                with pysam.VariantFile(str(folder / ('bench_' + name) / 'tp-base.vcf.gz')) as stream:
                    sets.append({variant_key(row) for row in stream})
            changes[name] = dict(gained_truth=sorted(sets[1]-sets[0]), lost_truth=sorted(sets[0]-sets[1]))
        save(directory / 'truth_deltas.json', changes)
    print(json.dumps({name: {k: v for k, v in m.items() if k != 'strata'} for name, m in metrics.items()}, indent=2))


if __name__ == '__main__':
    main()
