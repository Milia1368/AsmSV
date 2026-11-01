"""Run-level sequence/structure statistics; all scoring parameters are explicit."""
from bisect import bisect_left, bisect_right
from asmsv_next.constants import MATCH_REWARD
from dataclasses import dataclass
import math
from typing import Optional
from asmsv_next.alignment_tags import SourcePreset, TagError, synchronize
from asmsv_next.score_support import Interval, stable_id
QUERY = {0, 1, 4, 7, 8}
REFERENCE = {0, 2, 3, 7, 8}
ALIGNED = {0, 7, 8}

def normalize_cigar(cigar):
    result = []
    for op, length in cigar:
        if type(op) is not int or op not in range(9) or type(length) is not int or (length <= 0):
            raise ValueError('Invalid CIGAR operation')
        if result and result[-1][0] == op:
            result[-1] = (op, result[-1][1] + length)
        else:
            result.append((op, length))
    if not result:
        raise ValueError('Missing CIGAR')
    for i, (op, _) in enumerate(result):
        if op == 5 and i not in (0, len(result) - 1):
            raise ValueError('Internal hard clip in CIGAR')
        if op == 4 and (not (i == 0 or i == len(result) - 1 or (i == 1 and result[0][0] == 5) or (i == len(result) - 2 and result[-1][0] == 5))):
            raise ValueError('Internal soft clip in CIGAR')
    return tuple(result)

@dataclass(frozen=True)
class RawAlignment:
    q_name: str
    ref_name: str
    ref_start: int
    ref_end: int
    cigar: tuple
    raw_cigar: tuple
    mapq: int
    flag: int
    strand: str
    qs: int
    qe: int
    query_length: int
    raw_as: Optional[float]
    tags: tuple
    source: SourcePreset
    source_id: str
    haplotype: str
    alignment_id: str
    record_id: str
    issues: tuple = ()

    @classmethod
    def from_dict(cls, data):
        raw_cigar = tuple((tuple(x) for x in data.get('cigar', ())))
        cigar = normalize_cigar(raw_cigar)
        left = 0
        for op, length in cigar:
            if op in (4, 5):
                left += length
            else:
                break
        total = sum((n for op, n in cigar if op in QUERY or op == 5))
        aligned = sum((n for op, n in cigar if op in ALIGNED or op == 1))
        strand = data.get('strand', '-' if data.get('is_reverse', False) else '+')
        if strand not in ('+', '-'):
            raise ValueError('Invalid alignment strand')
        qs = data.get('qs', left if strand == '+' else total - left - aligned)
        qe = data.get('qe', left + aligned if strand == '+' else total - left)
        if type(qs) is not int or type(qe) is not int or qs < 0 or (qe - qs != aligned):
            raise ValueError('CIGAR and query coordinates disagree')
        query_length = data.get('query_length', max(total, qe))
        if type(query_length) is not int or query_length < max(total, qe):
            raise ValueError('Invalid full query length')
        start = data.get('ref_start', 0)
        if type(start) is not int or start < 0:
            raise ValueError('Invalid reference start')
        end = start + sum((n for op, n in cigar if op in REFERENCE))
        if data.get('ref_end', end) != end or end <= start or aligned <= 0:
            raise ValueError('CIGAR and reference coordinates disagree')
        tags = dict(data.get('tags', {}))
        native = data.get('raw_as', tags.get('AS', data.get('as_score')))
        issues = []
        if native is not None and (type(native) not in (int, float) or not math.isfinite(native)):
            issues.append('INVALID_NATIVE_AS')
            native = None
        mapq = data.get('mapq', 0)
        if type(mapq) is not int or not 0 <= mapq <= 255:
            raise ValueError('Invalid MAPQ')
        source_data = data.get('source', {})
        if isinstance(source_data, SourcePreset):
            source = source_data
        else:
            source_data = dict(source_data)
            if 'overrides' in source_data:
                source_data['overrides'] = tuple(source_data['overrides'])
            source = SourcePreset(**source_data)
        name, chrom = (data.get('q_name', ''), data.get('ref_name', ''))
        if not name or not chrom:
            raise ValueError('Missing query/reference name')
        source_id = data.get('source_id', 'unspecified')
        haplotype = data.get('haplotype', 'unspecified')
        identity = stable_id('ALN', source_id, haplotype, name, chrom, start, end, strand, qs, qe, cigar)
        record_id = stable_id('RAW', identity, mapq, int(data.get('flag', 0)), [(k, repr(v)) for k, v in sorted(tags.items())], native)
        return cls(name, chrom, start, end, cigar, raw_cigar, mapq, int(data.get('flag', 0)), strand, qs, qe, query_length, native, tuple(sorted(tags.items())), source, source_id, haplotype, identity, record_id, tuple(issues))

@dataclass(frozen=True)
class Run:
    qs: int
    qe: int
    rs: int
    re: int
    op: int
    mismatch_fraction: Optional[float] = 0.0
    positional: bool = True
    origin: str = 'CIGAR'

@dataclass(frozen=True)
class GapEvent:
    event_id: str
    alignment_id: str
    op: int
    query_interval: Interval
    reference_interval: Interval
    cost: float
    structural: bool

@dataclass(frozen=True)
class ScoreSummary:
    sequence_score: Optional[float]
    evaluation_length: int
    internal_structure_cost: float
    matches: Optional[float]
    mismatches: Optional[float]
    positional: bool

    @property
    def aq(self):
        return None if self.sequence_score is None or not self.evaluation_length else self.sequence_score / self.evaluation_length

    @property
    def node_reward(self):
        return None if self.sequence_score is None else self.sequence_score - self.internal_structure_cost

def geometry(raw):
    left = 0
    for op, length in raw.cigar:
        if op in (4, 5):
            left += length
        else:
            break
    q = raw.qs - left if raw.strand == '+' else 0
    axis = raw.qe + left
    r, result = (raw.ref_start, [])
    for op, length in raw.cigar:
        q2 = q + (length if op in QUERY or op == 5 else 0)
        r2 = r + (length if op in REFERENCE else 0)
        s, e = (q, q2) if raw.strand == '+' else (axis - q2, axis - q)
        result.append(Run(s, e, r, r2, op))
        q, r = (q2, r2)
    return tuple(result)

def _part(run, offset, length, op, raw, fraction, positional, origin):
    rs, re = (run.rs + offset, run.rs + offset + length)
    qs, qe = (run.qs + offset, run.qs + offset + length) if raw.strand == '+' else (run.qe - offset - length, run.qe - offset)
    return Run(qs, qe, rs, re, op, fraction, positional, origin)

class RunIndex:
    """Prefix sums over positive-query runs plus point costs for deletions."""

    def __init__(self, runs, gaps, preset):
        ordered = tuple(runs)
        if all((a.qe <= b.qs for a, b in zip(ordered, ordered[1:]))):
            self.runs = ordered
        elif all((a.qs >= b.qe for a, b in zip(ordered, ordered[1:]))):
            self.runs = ordered[::-1]
        else:
            raise ValueError('Evidence runs overlap or are out of query order')
        self.starts = tuple((r.qs for r in self.runs))
        self.ends = tuple((r.qe for r in self.runs))
        self.gaps = tuple(gaps)
        self._query_gaps = tuple(sorted(enumerate(self.gaps), key=lambda x: (x[1].query_interval.start, x[0])))
        self._query_gap_starts = tuple((g.query_interval.start for _, g in self._query_gaps))
        self.preset = preset
        self.positive = tuple((r for r in self.runs if r.qe > r.qs and r.op in (0, 1, 7, 8)))
        self.positive_ends = tuple((r.qe for r in self.positive))
        self.coefficients, prefixes = ([], [[0.0] * 7])
        for run in self.positive:
            n = run.qe - run.qs
            if run.op in ALIGNED:
                f = run.mismatch_fraction
                row = [0.0 if f is None else MATCH_REWARD * (1 - f) - preset.mismatch * f, 0.0, 1.0, 0.0 if f is None else 1 - f, 0.0 if f is None else f, float(f is None), float(not run.positional)]
            elif n >= preset.structural_gap_min:
                row = [0.0, preset.gap_cost(n) / n, 0.0, 0.0, 0.0, 0.0, 0.0]
            else:
                row = [-preset.gap_cost(n) / n, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
            self.coefficients.append(tuple(row))
            prefixes.append([a + n * b for a, b in zip(prefixes[-1], row)])
        self.prefixes = tuple((tuple(row) for row in prefixes))
        deletions = [g for g in gaps if g.op in (2, 3)]
        if any((a.query_interval.start > b.query_interval.start for a, b in zip(deletions, deletions[1:]))):
            deletions.reverse()
        self.deletion_positions = tuple((g.query_interval.start for g in deletions))
        points = [(0.0, 0.0, 0.0)]
        for gap in deletions:
            length = gap.reference_interval.length
            delta = (0.0, gap.cost, 0.0) if gap.structural else (-gap.cost, 0.0, float(length))
            points.append(tuple((a + b for a, b in zip(points[-1], delta))))
        self.deletion_prefix = tuple(points)

    def cumulative(self, q):
        i = bisect_right(self.positive_ends, q)
        result = list(self.prefixes[i])
        if i < len(self.positive):
            run = self.positive[i]
            length = max(0, min(q, run.qe) - run.qs)
            result = [a + length * b for a, b in zip(result, self.coefficients[i])]
        return result

    def score(self, start, end):
        if start > end:
            raise ValueError('Reversed score interval')
        row = [b - a for a, b in zip(self.cumulative(start), self.cumulative(end))]
        lo = bisect_right(self.deletion_positions, start)
        hi = bisect_left(self.deletion_positions, end)
        if lo < hi:
            delta = [b - a for a, b in zip(self.deletion_prefix[lo], self.deletion_prefix[hi])]
            row[0] += delta[0]
            row[1] += delta[1]
            row[2] += delta[2]
        known = row[5] < 0.5
        return ScoreSummary(row[0] if known else None, int(round(row[2])), row[1], row[3] if known else None, row[4] if known else None, known and row[6] < 0.5)

    def local_runs(self, start, end):
        lo = max(0, bisect_right(self.ends, start) - 1)
        hi = bisect_right(self.starts, end)
        return self.runs[lo:hi]

    def local_gaps(self, start, end):
        lo = bisect_left(self._query_gap_starts, start)
        hi = bisect_right(self._query_gap_starts, end)
        return tuple((g for _, g in sorted(self._query_gaps[lo:hi], key=lambda x: x[0])))

class EvidenceRoot:

    def __init__(self, raw, preset):
        self.raw, self.preset = (raw, preset)
        self.geometry = geometry(raw)
        indexed = tuple(enumerate(self.geometry))
        self.geometry_query_order = indexed if raw.strand == '+' else indexed[::-1]
        self.geometry_starts = tuple((run.qs for _, run in self.geometry_query_order))
        self.geometry_ends = tuple((run.qe for _, run in self.geometry_query_order))
        tags, issues = (dict(raw.tags), list(raw.issues))
        annotations, origin = (None, 'EQX')
        for kind in ('MD', 'cs'):
            if kind not in tags:
                continue
            try:
                parsed = synchronize(raw.cigar, tags[kind], kind)
                canonical = []
                for chunks in parsed:
                    merged = []
                    for symbol, n in chunks or ():
                        if merged and merged[-1][0] == symbol:
                            merged[-1] = (symbol, merged[-1][1] + n)
                        else:
                            merged.append((symbol, n))
                    canonical.append(merged if chunks is not None else None)
                if annotations is None:
                    annotations, origin = (canonical, kind)
                elif canonical != annotations:
                    issues.append('MD_CS_CONFLICT')
            except (TagError, TypeError, ValueError):
                issues.append('INVALID_' + kind.upper())
        m_bases = sum((n for op, n in raw.cigar if op == 0))
        explicit_x = sum((n for op, n in raw.cigar if op == 8))
        indels = sum((n for op, n in raw.cigar if op in (1, 2)))
        nm = tags.get('NM')
        nm_valid = type(nm) is int and nm >= 0 and (0 <= nm - indels - explicit_x <= m_bases)
        if 'NM' in tags and (not nm_valid):
            issues.append('NM_UNRECONCILED' if type(nm) is int and nm >= 0 else 'INVALID_NM')
        fraction = (nm - indels - explicit_x) / m_bases if nm_valid and m_bases else None
        runs, gaps = ([], [])
        for index, run in enumerate(self.geometry):
            op, length = raw.cigar[index]
            if op in ALIGNED:
                if annotations is not None:
                    offset = 0
                    for symbol, n in annotations[index]:
                        runs.append(_part(run, offset, n, 7 if symbol == '=' else 8, raw, float(symbol == 'X'), True, origin))
                        offset += n
                else:
                    runs.append(Run(run.qs, run.qe, run.rs, run.re, op, 0.0 if op == 7 else 1.0 if op == 8 else fraction, op != 0, 'EQX' if op != 0 else 'NM' if nm_valid else 'UNKNOWN'))
            elif op in (1, 2, 3):
                runs.append(run)
                if op == 3:
                    issues.append('REFERENCE_SKIP')
                gaps.append(GapEvent(stable_id('GAP', raw.alignment_id, index), raw.alignment_id, op, Interval(run.qs, run.qe), Interval(run.rs, run.re), preset.gap_cost(length), length >= preset.structural_gap_min))
        self.index = RunIndex(runs, gaps, preset)
        summary = self.index.score(raw.qs, raw.qe)
        if summary.positional and nm_valid:
            observed = sum((r.qe - r.qs for r in runs if r.op == 8)) + indels
            if observed != nm:
                issues.append('NM_UNRECONCILED')
        self.issues = tuple(sorted(set(issues)))

    def mapped_bounds(self, start, end):
        runs = [r for r in self.index.local_runs(start, end) if r.op in ALIGNED and min(end, r.qe) > max(start, r.qs)]
        if not runs:
            return None
        return (max(start, min((r.qs for r in runs))), min(end, max((r.qe for r in runs))))

    def cigar_view(self, start, end):
        raw, cigar, spans = (self.raw, [], [])
        lo = max(0, bisect_right(self.geometry_ends, start) - 1)
        hi = bisect_right(self.geometry_starts, end)
        local = self.geometry_query_order[lo:hi]
        if raw.strand == '-':
            local = local[::-1]
        for _, run in local:
            if run.op in ALIGNED or run.op == 1:
                s, e = (max(start, run.qs), min(end, run.qe))
                if s >= e:
                    continue
                cigar.append((run.op, e - s))
                if run.op in ALIGNED:
                    a, b = (run.rs + s - run.qs, run.rs + e - run.qs) if raw.strand == '+' else (run.re - (e - run.qs), run.re - (s - run.qs))
                    spans.append((a, b))
            elif run.op in (2, 3) and start < run.qs < end:
                cigar.append((run.op, run.re - run.rs))
        if not spans:
            raise ValueError('View has no mapped sequence')
        before = start if raw.strand == '+' else raw.query_length - end
        after = raw.query_length - end if raw.strand == '+' else start
        if before:
            cigar.insert(0, (5, before))
        if after:
            cigar.append((5, after))
        return (list(normalize_cigar(cigar)), min((a for a, b in spans)), max((b for a, b in spans)))
