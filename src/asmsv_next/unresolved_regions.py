"""Observational annotations for incomplete selected-path interpretations.

One region belongs to one detector context and one physical haplotype contig.
This module neither creates variant events nor changes atom ownership. Reference
coordinates alone never join regions or establish a mutation history.
"""
import hashlib
from collections import defaultdict

from asmsv_next.config import SVConfig

SCHEMA = 'asmsv-regions-v1'
POLICY = 'asmsv-unresolved-regions-v1'


def reason_categories(reasons):
    categories = set()
    for reason in reasons:
        if reason == 'ASSEMBLY_GAP':
            categories.add('ASSEMBLY_GAP')
        elif reason in ('REFERENCE_GAP', 'REFERENCE_UNAVAILABLE'):
            categories.add('REFERENCE_GAP')
        elif 'BUDGET' in reason:
            categories.add('COMPUTE_LIMIT')
        elif reason.startswith('OPEN_'):
            categories.add('OPEN_BOUNDARY')
        elif 'SEED_' in reason:
            categories.add('SEQUENCE_AMBIGUITY')
        elif any(word in reason for word in ('OWNERSHIP', 'UNMAPPED', 'OVERLAP', 'ANCHORS')):
            categories.add('ALIGNMENT_AMBIGUITY')
        elif 'RETENTION' in reason or reason == 'SOURCE_GROUPING_AMBIGUOUS':
            categories.add('SOURCE_UNCERTAINTY')
        else:
            categories.add('STRUCTURE_UNRESOLVED')
    return sorted(categories or {'STRUCTURE_UNRESOLVED'})


def _exact_anchor(block, point, outgoing, k):
    from asmsv_next.structural_detection import mapped_parts, map_port
    start, end = (point-k, point) if outgoing else (point, point+k)
    if not block.qs <= start < end <= block.qe:
        return False
    for qs, qe, rs, re in mapped_parts(block, (start, end), exact=True):
        if (qs, qe) != (start, end):
            continue
        expected = (re if block.strand == '+' else rs) if outgoing else (rs if block.strand == '+' else re)
        # One contiguous exact run is required. A zero-query deletion inside
        # the window or a two-valued boundary projection is not an exact port.
        return map_port(block, point) == expected
    return False


def _boundary(detector, index, outgoing):
    from asmsv_next.structural_detection import map_port
    if index is None:
        return dict(state='OPEN', block_index=None, query_position=None, port=None, exact_anchor=False)
    block = detector.chain[index]
    point = block.qe if outgoing else block.qs
    position = map_port(block, point)
    side = ('R' if block.strand == '+' else 'L') if outgoing else ('L' if block.strand == '+' else 'R')
    return dict(state='MAPPED' if position is not None else 'UNMAPPED', block_index=index,
                query_position=point, port=None if position is None else
                dict(chrom=block.ref_name, position=position, side=side),
                exact_anchor=_exact_anchor(block, point, outgoing, detector.k))


def _path_evidence(detector, context):
    """Retain certified local joins; record every rejected transition explicitly."""
    from asmsv_next.structural_detection import identifier, junction, map_port
    left, right = context['left_backbone'], context['right_backbone']
    indices = ([left] if left is not None else []) + list(context['member_block_indices'])
    if right is not None and right not in indices:
        indices.append(right)
    bounds = {i: [detector.chain[i].qs, detector.chain[i].qe] for i in indices}
    conflicts = set()
    for i, j in zip(indices, indices[1:]):
        if bounds[j][0] < bounds[i][1]:
            if detector.exact_overlap(detector.chain[i], detector.chain[j], bounds[j][0], bounds[i][1]):
                bounds[j][0] = bounds[i][1]
            else:
                conflicts.add((i, j))
    occurrences = []
    for index in indices:
        block = detector.chain[index]
        qs, qe = bounds[index]
        a, b = map_port(block, qs), map_port(block, qe)
        occurrences.append(dict(occurrence_id=identifier('OCC', detector.haplotype, detector.name,
            block.score_root.raw.alignment_id, qs, qe), block_index=index,
            alignment_id=block.score_root.raw.alignment_id, query=[qs, qe],
            reference=[block.ref_name, min(a, b) if a is not None and b is not None else block.ref_start,
                       max(a, b) if a is not None and b is not None else block.ref_end, block.strand],
            role='LEFT_ANCHOR' if index == left else 'RIGHT_ANCHOR' if index == right else 'MEMBER',
            entry=None if a is None else dict(chrom=block.ref_name, position=a, side='L' if block.strand == '+' else 'R'),
            exit=None if b is None else dict(chrom=block.ref_name, position=b, side='R' if block.strand == '+' else 'L')))
    edges, breaks = [], []
    for n, (i, j) in enumerate(zip(indices, indices[1:])):
        a, b = detector.chain[i], detector.chain[j]
        qs, qe = bounds[i][1], bounds[j][0]
        first, second = occurrences[n], occurrences[n+1]
        reason = None
        if (i, j) in conflicts or qe < qs or bounds[i][0] >= qs or qe >= bounds[j][1]:
            reason = 'QUERY_OWNERSHIP_CONFLICT'
        elif first['exit'] is None or second['entry'] is None:
            reason = 'UNMAPPED_PORT'
        elif not (_exact_anchor(a, qs, True, detector.k) and _exact_anchor(b, qe, False, detector.k)):
            reason = 'INSUFFICIENT_EXACT_ANCHOR'
        elif qe-qs > SVConfig.SEED_MAX_BASES:
            reason = 'JUNCTION_SEQUENCE_BUDGET'
        elif any(base not in 'ACGT' for base in detector.sequence[qs:qe].upper()):
            reason = 'ASSEMBLY_GAP'
        same_axis = a.ref_name == b.ref_name and a.strand == b.strand
        rg = (1 if a.strand == '+' else -1) * (second['entry']['position']-first['exit']['position']) if same_axis and first['exit'] and second['entry'] else None
        if reason is None and same_axis:
            if min(qe-qs, rg) >= SVConfig.MIN_SV_SIZE:
                reason = 'TWO_SIDED_RESIDUAL'
            elif abs((qe-qs)-rg) < SVConfig.MIN_SV_SIZE:
                # A normal selected-path step is not a novel adjacency.
                continue
        if reason:
            breaks.append(dict(left_occurrence=first['occurrence_id'], right_occurrence=second['occurrence_id'],
                               query=[qs, qe], reason=reason))
            continue
        edge = junction(detector.haplotype, detector.name, context['context_id'], first['exit'], second['entry'], [qs, qe])
        edge['junction_id'] = identifier('JUNC', context['context_id'], first['occurrence_id'], second['occurrence_id'], edge['left'], edge['right'])
        inserted = detector.sequence[qs:qe].upper()
        edge.update(evidence='SELECTED_PATH_EXACT_PORTS', adjacency_status='OBSERVED',
                    left_occurrence=first['occurrence_id'], right_occurrence=second['occurrence_id'],
                    inserted_sequence=inserted, sequence_length=len(inserted),
                    sequence_sha256=hashlib.sha256(inserted.encode()).hexdigest(), noncollinear=not same_axis or rg < 0)
        edges.append(edge)
    return occurrences, edges, breaks


def annotate_regions(detector, calls):
    from asmsv_next.structural_detection import identifier
    if not any(context['resolution'] != 'RESOLVED' for context in detector.contexts):
        return []
    regions = []
    calls_by_context, calls_by_id = defaultdict(list), defaultdict(list)
    events_by_context, events_by_atom = defaultdict(set), defaultdict(set)
    event_index = {event['event_id']: event for event in detector.events}
    for candidate in calls:
        for context_id in str(candidate.info.get('CONTEXT_ID', '')).split(','):
            calls_by_context[context_id].append(candidate)
        if candidate.id:
            calls_by_id[candidate.id].append(candidate)
    for event in detector.events:
        events_by_context[event.get('context_id')].add(event['event_id'])
        for atom_id in event.get('components', []) + ([event['basis_atom']] if event.get('basis_atom') else []):
            events_by_atom[atom_id].add(event['event_id'])
    for context in detector.contexts:
        if context['resolution'] == 'RESOLVED':
            continue
        region_id = identifier('REGION', POLICY, context['context_id'])
        context['region_id'] = region_id
        occurrences, edges, breaks = _path_evidence(detector, context)
        # Existing template/atom edges remain available. Add each physical join
        # once; never let an uncertain transition replace a certified edge.
        existing = context.setdefault('junctions', [])
        def geometry(edge):
            return (tuple(sorted((p['chrom'], p['position'], p['side']) for p in (edge['left'], edge['right']))),
                    tuple(edge.get('query_support') or []))
        by_geometry = {geometry(edge): edge for edge in existing}
        for edge in edges:
            previous = by_geometry.get(geometry(edge))
            if previous is None:
                existing.append(edge)
                by_geometry[geometry(edge)] = edge
            else:
                previous.update({key: value for key, value in edge.items() if key != 'junction_id'})
        if existing and context['resolution'] == 'UNRESOLVED':
            context['resolution'] = 'PARTIAL'
        for edge in existing:
            edge.setdefault('evidence', 'TEMPLATE_RECORDED')
            edge.setdefault('adjacency_status', 'TEMPLATE_RECORDED')
        boundaries = dict(left=_boundary(detector, context['left_backbone'], True),
                          right=_boundary(detector, context['right_backbone'], False))
        qleft = boundaries['left']['query_position']
        qright = boundaries['right']['query_position']
        qleft = 0 if qleft is None else qleft
        qright = len(detector.sequence) if qright is None else qright
        qcoords = [qleft, qright] + [p for o in occurrences if o['role'] == 'MEMBER' for p in o['query']]
        qcoords += context.get('query_span', [])
        qspan = [min(qcoords), max(qcoords)]
        intervals = []
        def interval(chrom, start, end, role, strand='.'):
            item = dict(chrom=chrom, start=start, end=end, role=role, strand=strand)
            if item not in intervals:
                intervals.append(item)
        lp, rp = boundaries['left']['port'], boundaries['right']['port']
        if lp and rp and lp['chrom'] == rp['chrom']:
            interval(lp['chrom'], min(lp['position'], rp['position']), max(lp['position'], rp['position']), 'TARGET')
        else:
            for port in (lp, rp):
                if port:
                    interval(port['chrom'], port['position'], port['position'], 'BOUNDARY_PORT')
        for occurrence in occurrences:
            if occurrence['role'] == 'MEMBER':
                chrom, start, end, strand = occurrence['reference']
                interval(chrom, start, end, 'MEMBER', strand)
        for chrom, start, end, strand in context.get('source_intervals', []):
            interval(chrom, start, end, 'SOURCE', strand)
        atom_ids = set(context['atom_ids'])
        indices = {o['block_index'] for o in occurrences}
        for index in indices:
            for atom_id in detector.cache.get(index, []):
                atom = detector.atoms[atom_id]
                if qspan[0] <= atom.query[0] <= atom.query[1] <= qspan[1]:
                    atom_ids.add(atom_id)
        event_ids = set(events_by_context[context['context_id']])
        for atom_id in atom_ids:
            event_ids.update(events_by_atom[atom_id])
        event_ids = sorted(event_ids)
        # CRR is an optional observational flag: at least two certified
        # noncollinear joins AND at least two incompatible member source runs.
        # Fragmentation of one monotone inversion source is not sufficient.
        members = [o for o in occurrences if o['role'] == 'MEMBER']
        incompatible = False
        for a, b in zip(members, members[1:]):
            ac, ar0, ar1, astrand = a['reference']
            bc, br0, br1, bstrand = b['reference']
            if ac != bc or astrand != bstrand or (br1 <= ar0 if astrand == '+' else br0 >= ar1):
                incompatible = True
        crr = incompatible and sum(edge.get('noncollinear', False) for edge in edges) >= 2
        labels = ['UNRESOLVED_REGION'] + (['CRR'] if crr else [])
        context['region_annotations'] = labels
        region = dict(record_kind='UNRESOLVED_REGION', schema=SCHEMA, policy=POLICY,
            region_id=region_id, context_id=context['context_id'], query_name=detector.name,
            haplotype=detector.haplotype, resolution=context['resolution'], labels=labels,
            reasons=context['issues'] or ['INCOMPLETE_INTERPRETATION'],
            reason_categories=reason_categories(context['issues']), evidence_scope='SELECTED_PATH',
            alt_assessment='NOT_AVAILABLE', query_span=qspan, boundary_order='FORWARD' if qleft <= qright else 'OVERLAP',
            boundaries=boundaries, reference_intervals=intervals,
            members=dict(block_indices=list(context['member_block_indices']), atom_ids=sorted(atom_ids),
                         event_ids=event_ids, junction_ids=sorted(e['junction_id'] for e in existing)),
            observed_path=occurrences, junctions=existing, unretained_transitions=breaks,
            separators=context['separators'], sequence_checks=context.get('sequence_checks', []),
            counts_as_variant=False, genotype=None)
        regions.append(region)
        for atom_id in atom_ids:
            atom = detector.atoms.get(atom_id)
            if atom:
                _add_region(atom.candidate.info, region_id)
        for candidate in calls_by_context[context['context_id']]:
            _add_region(candidate.info, region_id)
        for event_id in event_ids:
            for candidate in calls_by_id[event_id]:
                _add_region(candidate.info, region_id)
            event_index[event_id].setdefault('region_ids', []).append(region_id)
    return regions


def _add_region(info, region_id):
    info['REGION_ID'] = ','.join(sorted(set(info.get('REGION_ID', '').split(',')) - {''} | {region_id}))
