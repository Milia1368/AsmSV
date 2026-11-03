"""Presentation-only REF/ALT schematics from already reported local structures.

These are large-edit cartoons, not base-level allele reconstruction. Incomplete
placement is explicit; no detector decision or benchmark endpoint is changed.
"""
from copy import deepcopy


def edit_context(atom):
    info = atom.get('info', {})
    return info.get('ALLELE_CONTEXT', info.get('EDIT_CONTEXT'))


def normalized_key(query, strand):
    return query[0] if strand == '+' else -query[1]


def edit_list(structure):
    result = []
    for atom in structure.get('atoms', []):
        kind = atom['type']
        if kind in ('DEL', 'DELINS'):
            result.append(dict(id=atom['atom_id'], type='DEL', chrom=atom['chrom'],
                               start=atom['pos'], end=atom['end'], length=atom['end']-atom['pos'],
                               query=atom['query'], context=edit_context(atom)))
        if kind in ('INS', 'DELINS'):
            result.append(dict(id=atom['atom_id'], type='INS', chrom=atom['chrom'],
                               start=atom['pos'], end=atom['pos'], length=atom['query'][1]-atom['query'][0],
                               query=atom['query'], context=edit_context(atom)))
    return result


def partition(chrom, start, end, edits):
    cuts = {start, end}
    for edit in edits:
        if edit['chrom'] == chrom:
            cuts.update(p for p in (edit['start'], edit['end']) if start < p < end)
    cuts = sorted(cuts)
    bins, sequence_number, deletion_number = [], 0, 0
    for left, right in zip(cuts, cuts[1:]):
        deleted = [e for e in edits if e['type']=='DEL' and e['chrom']==chrom
                   and e['start'] <= left and right <= e['end']]
        if deleted:
            deletion_number += 1
            label = 'D' if sum(e['type']=='DEL' for e in edits)==1 else 'D'+str(deletion_number)
        else:
            label = chr(65 + sequence_number) if sequence_number < 26 else 'S'+str(sequence_number+1)
            sequence_number += 1
        bins.append(dict(label=label, kind='segment', chrom=chrom, start=left, end=right,
                         length=right-left, direction='+', deletion_ids=[e['id'] for e in deleted]))
    if len(bins) == 1 and not bins[0]['deletion_ids']:
        bins[0]['label'] = 'S'
    return bins


def insertion(edit, label='X'):
    return dict(kind='insertion', label=label, length=edit['length'], direction=None,
                edit_id=edit['id'], context=edit['context'])


def altered_tokens(bins, edits, direction, query_strand='+'):
    tokens = deepcopy(bins if direction == '+' else list(reversed(bins)))
    deletions = {e['id'] for e in edits if e['type']=='DEL'}
    for token in tokens:
        token['direction'] = direction
        if deletions.intersection(token['deletion_ids']):
            token['kind'] = 'missing'
    insertions = sorted((e for e in edits if e['type']=='INS'),
                        key=lambda e:normalized_key(e['query'], query_strand))
    for number, edit in enumerate(insertions):
        cut = edit['start']
        boundary = next((i for i, b in enumerate(tokens) if b.get('kind') != 'insertion'
                         and (b['start'] if direction=='+' else b['end'])==cut), None)
        if boundary is None and tokens and cut == (bins[-1]['end'] if direction=='+' else bins[0]['start']):
            boundary = len(tokens)
        if boundary is None:
            raise ValueError('INSERTION_PLACEMENT_UNAVAILABLE')
        tokens.insert(boundary, insertion(edit, 'X' if len(insertions)==1 else 'X'+str(number+1)))
    return tokens


def group(role, title, tokens, direction='+', **extra):
    return dict(kind='group', role=role, title=title, tokens=tokens, direction=direction, **extra)


def diagram_for_structure(structure):
    event, context = structure['event'], structure['context']
    info = event.get('info', {})
    result = dict(status='partial', frame=info.get('FRAME', context.get('frame')),
                  haplotype=event.get('haplotype'), query_name=event.get('query_name'),
                  context_id=event.get('context_id'), summaries=[], projections=[],
                  ref=[], alt=[], notes=[])
    edits = edit_list(structure)
    names = {'REFERENCE':'发生在参考框架内', 'ADDED_DUPLICATE':'发生在新增拷贝内',
             'DUPLICATION_JUNCTION':'发生在复制连接处'}
    result['summaries'] = [dict(type=e['type'], length=e['length'],
                                context=e['context'], text=names.get(e['context'], '编辑归属见证据详情'))
                           for e in sorted(edits,key=lambda e:(e['type']!='DEL',e['chrom'],e['start']))]
    result['projections'] = context.get('insertion_projections', [])
    try:
        sources = context.get('source_intervals', [])
        strand = context.get('target_strand')
        if len(sources)!=1 or strand not in ('+', '-'):
            raise ValueError('SOURCE_OR_TARGET_ORIENTATION_UNAVAILABLE')
        chrom, start, end, _ = sources[0]
        frame = result['frame']
        if frame == 'DUP':
            if info.get('DUP_MODE')!='ADJACENT' or info.get('SOURCE_STATE')!='RETAINED':
                raise ValueError('ONLY_RETAINED_ADJACENT_COPY_LAYOUT_IS_SUPPORTED')
            target = info.get('TARGET_POS')
            if info.get('TARGET_CHR')!=chrom or target not in (start, end):
                raise ValueError('TARGET_IS_NOT_A_SOURCE_BOUNDARY')
            segments = event.get('segments', [])
            occurrence = info.get('DUP_OCCURRENCE')
            if len(segments)!=2 or occurrence not in (0, 1):
                raise ValueError('COPY_VISIT_IDENTITY_UNAVAILABLE')
            local = [e for e in edits if e['context']=='ADDED_DUPLICATE']
            junction = [e for e in edits if e['context']=='DUPLICATION_JUNCTION']
            if len(local)+len(junction)!=len(edits) or any(e['type']!='INS' for e in junction):
                raise ValueError('UNSUPPORTED_COPY_EDIT_CONTEXT')
            if any(e['chrom']!=chrom or not start <= e['start'] <= e['end'] <= end for e in local):
                raise ValueError('COPY_EDIT_OUTSIDE_SOURCE')
            bins = partition(chrom, start, end, local)
            result['ref'] = [group('reference', '参考来源', bins, source=[chrom,start,end])]
            ordered = sorted(enumerate(segments), key=lambda p:normalized_key(p[1]['query'], strand))
            added_index = next(i for i,(original,_) in enumerate(ordered) if original==occurrence)
            if added_index != (0 if target==start else 1):
                raise ValueError('COPY_ORDER_AND_TARGET_DISAGREE')
            visits = []
            for original, segment in ordered:
                role = 'added' if original==occurrence else 'retained'
                direction = '+' if segment['source'][3]==strand else '-'
                if role=='retained' and direction!='+':
                    raise ValueError('RETAINED_SOURCE_NOT_COLLINEAR')
                parts = altered_tokens(bins, local if role=='added' else [], direction, strand)
                visits.append((normalized_key(segment['query'],strand), group(role,
                    '新增拷贝' if role=='added' else '来源保留遍历', parts, direction,
                    query=segment['query'], measured_length=info.get('DUP_SEQUENCE_LEN') if role=='added' else None)))
            for number, edit in enumerate(junction):
                visits.append((normalized_key(edit['query'],strand), dict(insertion(edit,'X' if len(junction)==1 else 'X'+str(number+1)),
                    role='junction', title='复制连接处插入')))
            result['alt'] = [v for _,v in sorted(visits,key=lambda v:v[0])]
            if info.get('DUP_ORIENT') in ('DIRECT','INVERTED'):
                expected = '+' if info['DUP_ORIENT']=='DIRECT' else '-'
                if next(g for g in result['alt'] if g.get('role')=='added')['direction']!=expected:
                    raise ValueError('COPY_ORIENTATION_DISAGREES')
        elif frame == 'INV':
            if any(e['context']!='REFERENCE' for e in edits):
                raise ValueError('UNSUPPORTED_INVERSION_EDIT_CONTEXT')
            span = event['reference_span']
            if span[0]!=chrom or not span[1] <= start < end <= span[2]:
                raise ValueError('INVERSION_REFERENCE_BOUNDARIES_UNAVAILABLE')
            if any(e['chrom']!=chrom or not span[1] <= e['start'] <= e['end'] <= span[2] for e in edits):
                raise ValueError('EDIT_OUTSIDE_REPORTED_FRAME')
            # Every outside-INV reference base must be explained by a deletion.
            bins = partition(chrom, span[1], span[2], edits + [dict(type='CUT',chrom=chrom,start=start,end=end)])
            outside = [b for b in bins if b['end']<=start or b['start']>=end]
            if any(not b['deletion_ids'] for b in outside):
                raise ValueError('UNEXPLAINED_REFERENCE_FLANK')
            inner = [b for b in bins if start<=b['start'] and b['end']<=end]
            internal = [e for e in edits if start <= e['start'] <= e['end'] <= end]
            left = [b for b in outside if b['end']<=start]
            right = [b for b in outside if b['start']>=end]
            missing = lambda bs:[dict(b,kind='missing') for b in deepcopy(bs)]
            # Boundary insertions are placed by their actual query order.
            segments = event.get('segments', [])
            if not segments or any(s['source'][3]==strand for s in segments):
                raise ValueError('INVERSION_ORIENTATION_UNAVAILABLE')
            qmin=min(s['query'][0] for s in segments);qmax=max(s['query'][1] for s in segments)
            before, after, internal_insertions = [], [], []
            for e in [e for e in edits if e['type']=='INS']:
                if e['query'][1]<=qmin or e['query'][0]>=qmax:
                    preceding = e['query'][1]<=qmin if strand=='+' else e['query'][0]>=qmax
                    (before if preceding else after).append(e)
                elif start<=e['start']<=end:
                    internal_insertions.append(e)
                else:
                    raise ValueError('INSERTION_PLACEMENT_UNAVAILABLE')
            inside_edits=[e for e in internal if e['type']=='DEL']+internal_insertions
            alt = missing(left)+[insertion(e) for e in before]+altered_tokens(inner,inside_edits,'-',strand)
            alt += [insertion(e) for e in after]+missing(right)
            result['ref'] = [group('reference','参考框架',bins,source=span)]
            result['alt'] = [group('inverted','倒位及耦合编辑',alt,'-', measured_length=event['query_span'][1]-event['query_span'][0])]
        else:
            raise ValueError('FRAME_NOT_SUPPORTED_BY_SCHEMATIC')
        result['status']='supported'
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        result['reason']=str(exc)
        result['notes'].append('现有字段不足以可靠排布此结构；ALT 不自动补推，请查看下方原始证据。')
    return result


def native_copy_diagram(record, structure=None):
    """Source/target geometry is explicit even without two query visits in JSON."""
    info = dict(record['info'])
    context = structure.get('context', {}) if structure else {}
    event = structure.get('event', {}) if structure else {}
    if structure:
        info.update(event.get('info', {}))
    subtype = info.get('SUBTYPE', info.get('DUP_TYPE', ''))
    orientation = info.get('DUP_ORIENT')
    if not orientation and subtype in ('piDUP_FR', 'piDUP_RF'):
        orientation = 'INVERTED'
    result = dict(status='partial', frame='DUP', haplotype=event.get('haplotype', 'VCF'),
                  query_name=event.get('query_name'), context_id=event.get('context_id'),
                  ref=[], alt=[], summaries=[], projections=context.get('insertion_projections', []),
                  annotations=[], notes=[], native_type='DUP')
    try:
        chrom = info.get('SOURCE_CHR', record['chrom'])
        start, end = int(info.get('SOURCE_START', record['pos'])), int(info.get('SOURCE_END', record['end']))
        target_chrom, target = info['TARGET_CHR'], int(info['TARGET_POS'])
        if not start < end or orientation not in ('DIRECT', 'INVERTED'):
            raise ValueError('SOURCE_OR_COPY_ORIENTATION_UNAVAILABLE')
        direction = '+' if orientation=='DIRECT' else '-'
        bins = partition(chrom, start, end, [])
        reference = group('reference','参考来源', bins, source=[chrom,start,end])
        retained = group('retained','来源保留遍历',deepcopy(bins))
        added = group('added','新增拷贝',altered_tokens(bins,[],direction), direction)
        result.update(source=[chrom,start,end], target=[target_chrom,target], copy_direction=direction,
                      ref=[reference], annotations=[
                          '原始类型 DUP；' + ('反向复制' if direction=='-' else '同向复制') + '，来源保留，目标处新增拷贝。'])
        if chrom==target_chrom and target in (start,end):
            result['alt'] = [added,retained] if target==start else [retained,added]
        else:
            # Separate loci, not a fabricated continuous source-to-target path.
            result.update(layout='source_target', alt=[retained,added])
        result['status']='supported'
    except (KeyError,ValueError,TypeError) as exc:
        result['reason']=str(exc)
    return result


def replacement_diagram(record):
    info=record['info']
    result=dict(status='partial',frame='DELINS',native_type='DELINS',haplotype='VCF',
                query_name=None,context_id=None,ref=[],alt=[],summaries=[],projections=[],notes=[],
                annotations=['同一局部替换的 DEL＋INS 表达；不声称唯一编辑分解或两个独立突变。'])
    try:
        deleted,inserted=int(info['DEL_LEN']),int(info['INS_LEN'])
        if min(deleted,inserted)<=0 or record['end']-record['pos']!=deleted:
            raise ValueError('REPLACEMENT_LENGTHS_DISAGREE')
        token=dict(label='D',kind='segment',chrom=record['chrom'],start=record['pos'],end=record['end'],
                   length=deleted,direction='+',deletion_ids=['replacement'])
        edit=dict(id='replacement',length=inserted,context='REFERENCE')
        result.update(status='supported',ref=[group('reference','被替换的参考片段',[token],
                                                    source=[record['chrom'],record['pos'],record['end']])],
                      alt=[group('replacement','DEL＋INS · 局部替换',
                                 [dict(token,kind='missing'),insertion(edit)])],
                      summaries=[dict(type='DEL',length=deleted,text='局部替换中的参考缺失'),
                                 dict(type='INS',length=inserted,text='同一替换中的新增序列')])
    except (KeyError,ValueError,TypeError) as exc:
        result['reason']=str(exc)
    return result


def diagrams_for_event(event):
    if event['main']['type']=='DUP':
        return [native_copy_diagram(event['main'], st) for st in event.get('structures', [])] or [native_copy_diagram(event['main'])]
    if event['main']['type']=='DELINS':
        return [replacement_diagram(event['main'])]
    return [diagram_for_structure(st) for st in event.get('structures', [])]
