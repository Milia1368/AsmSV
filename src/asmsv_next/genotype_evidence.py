"""Reference-callability guards after (unchanged) allele clustering.

An assembly with no matching ALT is not necessarily a reference observation.
This module only qualifies genotypes; it never creates or removes SV calls.
"""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict

from asmsv_next.config import SVConfig


def contig_evidence(chain, backbone_indices, query):
    """Reuse the AQ run index, without parsing an alignment CIGAR again.

    A reference traversal stops at an SV-sized indel, reference skip, or block
    boundary. Subthreshold edits and substitutions retain SV-level callability.
    Coordinates here are 0-based, half-open.
    """
    reference = []
    phase = defaultdict(list)
    for i in sorted(backbone_indices):
        block = chain[i]
        phase[block.ref_name].append((block.ref_start, block.ref_end))
        spans = []
        start = end = None
        runs = block.score_root.index.local_runs(block.qs, block.qe)
        for run in (runs if block.strand=='+' else reversed(runs)):
            if run.op in (2,3) and not block.qs < run.qs < block.qe:
                continue
            if run.op == 1 and min(run.qe,block.qe) <= max(run.qs,block.qs):
                continue
            if run.op in (0, 7, 8):
                lo, hi = max(run.qs, block.qs), min(run.qe, block.qe)
                if hi <= lo:
                    continue
                rs = run.rs + (lo-run.qs if block.strand == '+' else run.qe-hi)
                re = rs + hi-lo
                if start is None:
                    start, end = rs, re
                else:
                    end = max(end, re)
            elif run.op == 3 or run.op in (1, 2) and max(run.qe-run.qs, run.re-run.rs) >= SVConfig.MIN_SV_SIZE:
                if start is not None:
                    spans.append((start, end))
                start = end = None
            elif run.op == 2 and start is not None:
                end = min(block.ref_end,max(end, run.re))
        if start is not None:
            spans.append((start, end))
        reference.extend([block.ref_name, s, e, query, block.score_root.raw.alignment_id]
                         for s, e in spans if e > s)
    return dict(reference=reference, phase=[[chrom, min(s for s, e in spans),
                max(e for s, e in spans), query] for chrom, spans in sorted(phase.items())])


class IntervalIndex:
    def __init__(self, rows):
        grouped = defaultdict(list)
        for row in rows:
            grouped[row[0]].append(tuple(row[1:]))
        self.data = {}
        for chrom, intervals in grouped.items():
            intervals.sort()
            maximum, prefix = -1, []
            for row in intervals:
                maximum = max(maximum, row[1]); prefix.append(maximum)
            self.data[chrom] = (intervals, [r[0] for r in intervals], prefix)

    def overlaps(self, chrom, start, end):
        rows, starts, prefix = self.data.get(chrom, ([], [], []))
        stop = bisect_left(starts, end)
        begin = bisect_right(prefix, start, 0, stop)
        return [r for r in rows[begin:stop] if r[1] > start]


def genotype_regions(candidate):
    """DUP dosage is assessed at its target, not by donor-locus coverage."""
    c = candidate
    if c.type == 'DUP' or c.type == 'CPX' and 'DUP' in str(c.info.get('CPX_COMPONENTS', '')).split(','):
        if 'TARGET_POS' in c.info:
            p = int(c.info['TARGET_POS'])
            return [(c.info.get('TARGET_CHR', c.chrom), p, p+1)]
    if c.type in ('BND', 'TRA'):
        p = c.pos
        regions = [(c.chrom, max(0, p-1), p)]
        if c.info.get('CHR2') and c.info.get('POS2'):
            p2 = int(c.info['POS2'])
            regions.append((c.info['CHR2'], max(0, p2-1), p2))
        return sorted(regions)
    return [(c.chrom, c.pos, max(c.pos+1, c.end))]


class GenotypeEvidence:
    def __init__(self, bundles, flank=16):
        self.flank = flank
        reference, phase = [[], []], [[], []]
        for bundle in bundles:
            hap = {'hap1': 0, 'hap2': 1}[bundle['haplotype']]
            data = bundle.get('genotyping', {})
            reference[hap].extend(data.get('reference', []))
            phase[hap].extend(data.get('phase', []))
        self.reference = [IntervalIndex(rows) for rows in reference]
        self.phase = [IntervalIndex(rows) for rows in phase]
        # PS is chromosome-local, unique for each physical contig pair. Two
        # disconnected alignment spans still share phase within that pair.
        pairs = defaultdict(dict)
        for chrom, start, end, name in sorted(phase[0]):
            for s2, e2, name2 in self.phase[1].overlaps(chrom, start, end):
                key = (name, name2)
                anchor = max(start, s2)+1
                pairs[chrom][key] = min(anchor, pairs[chrom].get(key, anchor))
        self.ps = {}
        for chrom, anchors in sorted(pairs.items()):
            used = set()
            for pair, anchor in sorted(anchors.items(), key=lambda item: (item[1], item[0])):
                while anchor in used:
                    anchor += 1
                used.add(anchor); self.ps[(chrom, *pair)] = anchor

    def reference_state(self, hap, regions):
        names = set()
        for chrom, start, end in regions:
            start, end = max(0, start-self.flank), end+self.flank
            hits = self.reference[hap].overlaps(chrom, start, end)
            covering = [r for r in hits if r[0] <= start and r[1] >= end]
            # Partial competing alignments also preclude a unique traversal.
            identities = {(r[2], r[3]) for r in hits}
            if len(identities) > 1:
                return None, 'AMBIGUOUS_COVERAGE'
            if not covering:
                return None, 'NO_REFERENCE_SUPPORT'
            names.update(r[2] for r in covering)
        return names, 'PASS'

    def apply(self, calls, policy='evidence'):
        if policy not in ('assembly','evidence'):
            raise ValueError('Unknown genotype evidence policy: '+policy)
        for c in calls:
            alt = ['Hap1' in c.sources, 'Hap2' in c.sources]
            c.format_data = {}
            if not any(alt):
                c.genotype = './.'; c.format_data['FT'] = 'NO_ALT_SUPPORT'
                continue
            if all(alt):
                c.genotype = '1/1'
                continue
            other = 1 if alt[0] else 0
            regions = genotype_regions(c)
            copy_local = c.info.get('ROLE') == 'COMPONENT' and c.info.get('ALLELE_CONTEXT') in ('ADDED_DUPLICATE', 'DUPLICATION_JUNCTION')
            if copy_local:
                names, reason = None, 'NONREFERENCE_COMPONENT'
            else:
                names, reason = self.reference_state(other, regions)
            if names is None:
                # Preserve the established assembly dosage estimate by default,
                # but explicitly filter its unsupported reference assumption.
                # Strict evidence mode reports a missing counterpart instead.
                c.genotype = '1/.' if policy=='evidence' or copy_local else '0/1'
                c.format_data['FT'] = reason
                continue
            alleles = ['1' if present else '0' for present in alt]
            # A phase label belongs to the VCF chromosome. Cross-chromosome
            # copy targets and multi-port junctions retain dosage without phase.
            phase_names = [set(), set()]
            supported = 0 if alt[0] else 1
            if len(regions) == 1 and regions[0][0] == c.chrom:
                chrom, start, end = regions[0]
                hits = self.phase[supported].overlaps(chrom, start, end)
                phase_names[supported] = {r[2] for r in hits if r[2] in c.read_names}
                phase_names[other] = names
            if all(len(n) == 1 for n in phase_names):
                key = (c.chrom, next(iter(phase_names[0])), next(iter(phase_names[1])))
                ps = self.ps.get(key)
                if ps is not None:
                    c.genotype = '|'.join(alleles); c.format_data['PS'] = ps
                    continue
            c.genotype = '/'.join(sorted(alleles))
        return calls


def summarize_genotypes(calls):
    types = defaultdict(Counter)
    genotypes, reasons = Counter(), Counter()
    for c in calls:
        genotypes[c.genotype] += 1; types[c.type][c.genotype] += 1
        if c.format_data.get('FT'):
            reasons[c.format_data['FT']] += 1
    total = len(calls)
    partial = sum(n for gt, n in genotypes.items() if '.' in gt)
    filtered=sum(reasons.values())
    return dict(records=total, fully_called=total-partial, partial_or_missing=partial,
                call_rate=(total-partial)/total if total else None,
                genotype_filtered=filtered,unfiltered_records=total-filtered,
                unfiltered_fraction=(total-filtered)/total if total else None,
                genotypes=dict(sorted(genotypes.items())), by_type={k: dict(sorted(v.items())) for k, v in sorted(types.items())},
                filter_reasons=dict(sorted(reasons.items())))
