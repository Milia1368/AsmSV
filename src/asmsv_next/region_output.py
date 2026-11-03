"""Standalone unresolved-region JSONL and BED projections, without SV counting."""
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
from urllib.parse import quote
import uuid

from asmsv_next.unresolved_regions import POLICY, SCHEMA


def _regions(structure_dir, region_source):
    if region_source is not None:
        with Path(region_source).open() as stream:
            for line in stream:
                yield json.loads(line)
    else:
        for path in sorted(Path(structure_dir).glob('*.json')):
            yield from json.loads(path.read_text()).get('regions', [])


def write_regions(structure_dir, record_maps, output, aliases, contig_lengths, mask_applied=False, record_links=None, region_source=None):
    """Write one region per incomplete context; do not merge haplotypes or loci.

    Exact (view, record_id) references describe only records actually written.
    An empty list is valid, including when all member calls were BED-masked.
    """
    links = {key:list(value) for key,value in record_links.items()} if record_links is not None else defaultdict(list)
    if record_links is None:
        for path in record_maps:
            with Path(path).open() as stream:
                for line in stream:
                    record = json.loads(line)
                    for region_id in str(record['info'].get('REGION_ID', '')).split(','):
                        if region_id:
                            links[region_id].append(dict(view=record['view'], record_id=record['record_id']))
    target = Path(str(output)+'.regions.jsonl')
    bed = Path(str(output)+'.regions.bed')
    token = uuid.uuid4().hex
    temporary = [p.with_name(p.name+'.'+token+'.tmp') for p in (target, bed)]
    summary = Counter(regions=0, PARTIAL=0, UNRESOLVED=0, CRR=0, junctions=0)
    reasons = Counter()
    bed_rows = []
    metadata = dict(record_kind='REGION_METADATA', schema=SCHEMA, policy=POLICY,
        coordinate_system='0-based-half-open', evidence_scope='SELECTED_PATH',
        grouping='ONE_INCOMPLETE_CONTEXT_PER_HAPLOTYPE_CONTIG', counts_as_variant=False,
        candidate_mask_applied=mask_applied, region_mask_policy='RETAIN_ALL_DIAGNOSTICS',
        bed_port_projection='[cut,cut+1), or [length-1,length) at chromosome end',
        vcf_record_links='Exact view and record_id pairs; no record is implied when empty')
    try:
        with temporary[0].open('w') as stream:
            stream.write(json.dumps(metadata, sort_keys=True)+'\n')
            for original in _regions(structure_dir, region_source):
                region = dict(original)
                references = links.pop(region['region_id'], [])
                region['vcf_records'] = sorted((dict(view=item[0], record_id=item[1]) if isinstance(item, tuple) else item
                    for item in references), key=lambda x: (x['view'], x['record_id']))
                region['canonical_event_ids'] = sorted({aliases[key] for key in region['members']['event_ids'] if key in aliases})
                summary['regions'] += 1
                summary[region['resolution']] += 1
                summary['CRR'] += int('CRR' in region['labels'])
                summary['junctions'] += len(region['junctions'])
                reasons.update(region['reasons'])
                for interval in region['reference_intervals']:
                    chrom, start, end = interval['chrom'], interval['start'], interval['end']
                    length = contig_lengths[chrom]
                    if not 0 <= start <= end <= length:
                        raise ValueError(f"{region['region_id']}: invalid reference interval {interval}")
                    if start == end:
                        if length == 0:
                            continue
                        start = min(start, length-1)
                        end = start+1
                    values = (chrom, start, end, region['region_id'], 0, interval['strand'],
                              interval['role'], region['haplotype'], quote(region['query_name'], safe=''),
                              region['resolution'], ','.join(region['labels']), ','.join(region['reasons']))
                    bed_rows.append(values)
                stream.write(json.dumps(region, ensure_ascii=False, sort_keys=True, allow_nan=False)+'\n')
            if links:
                raise ValueError('VCF records reference absent unresolved regions: '+','.join(sorted(links)))
        with temporary[1].open('w') as stream:
            stream.write('#chrom\tstart\tend\tregion_id\tscore\tstrand\tinterval_role\thaplotype\tquery_name\tresolution\tlabels\treasons\n')
            for row in sorted(set(bed_rows), key=lambda row: (row[0], row[1], row[2], row[3], row[6])):
                stream.write('\t'.join(map(str, row))+'\n')
        for source, destination in zip(temporary, (target, bed)):
            os.replace(source, destination)
    finally:
        for path in temporary:
            if path.exists():
                path.unlink()
    return dict(jsonl=str(target), bed=str(bed), counts=dict(summary), reasons=dict(sorted(reasons.items())))
