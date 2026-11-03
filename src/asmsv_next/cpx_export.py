"""Read-only composite SV summaries: VERIX input and an offline HTML browser.

No detection, regrouping, sequence reconstruction or genotype imputation occurs.
Only Python's standard library is needed by this module.
"""
import argparse
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import tempfile
from urllib.parse import quote

from asmsv_next.cpx_diagram import diagrams_for_event
from asmsv_next.csv_scope import csv_class, composite_class, POLICY


LOG = logging.getLogger('asmSV-cpx')
SCHEMA = 'asmsv-csv-view-v1'
FIELDS = {
    'ASMSV_CSV_CLASS': ('1', 'String', 'Export selection stratum; original SVTYPE is retained; asmsv-csv-selection-v1'),
    'BKPS': ('.', 'String', 'VERIX single format TYPE:length-chrom:cut[-chrom:cut]; cuts are 0-based interbase; not an ALT path'),
    'ASMSV_COMPONENTS': ('.', 'String', 'id|type|chrom|VCF_POS|END|SVLEN|EDIT_CONTEXT; subfields percent-encoded'),
    'ASMSV_COMPONENT_GT': ('.', 'String', 'id|sample|GT|PS|FT; original component values; subfields percent-encoded'),
    'ASMSV_EXPORT_WARN': ('.', 'String', 'Limitations of this export; not caller FILTER values'),
}


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def parse_info(text):
    result = {}
    for token in text.split(';'):
        if token and token != '.':
            key, sep, value = token.partition('=')
            result[key] = value if sep else True
    return result


def values(value):
    return [] if value in (None, '.', '') else str(value).split(',')


def records(path):
    with open_text(path) as stream:
        columns = None
        for number, line in enumerate(stream, 1):
            if line.startswith('#CHROM\t'):
                columns = line.rstrip('\r\n').split('\t')
            elif line.startswith('#'):
                continue
            elif line.strip():
                fields = line.rstrip('\r\n').split('\t')
                if columns is None or len(fields) != len(columns) or len(fields) < 8:
                    raise ValueError(f'{path}:{number}: invalid VCF row/header')
                yield fields, parse_info(fields[7]), columns
        if columns is None:
            raise ValueError(f'{path}: missing #CHROM header')


def samples(fields, columns):
    keys = fields[8].split(':') if len(fields) > 8 else []
    return {name: dict(zip(keys, data.split(':'))) for name, data in zip(columns[9:], fields[9:])}


def compact_record(fields, info, columns):
    return dict(id=fields[2], chrom=fields[0], pos=int(fields[1]),
                end=int(info.get('END', int(fields[1]) + len(fields[3]) - 1)),
                type=info.get('SVTYPE', 'NA'), filter=fields[6],
                info={k: v for k, v in info.items() if k != 'SEQ'},
                samples=samples(fields, columns), ref_length=len(fields[3]),
                alt_lengths=[len(a) for a in fields[4].split(',')],
                shared_parents=values(info.get('PARENT_EVENT')))


def companion(vcf, suffix):
    candidates = [Path(str(vcf) + suffix), Path(str(vcf) + suffix + '.gz')]
    if str(vcf).endswith('.gz'):
        candidates.extend([Path(str(vcf)[:-3] + suffix), Path(str(vcf)[:-3] + suffix + '.gz')])
    return next((p for p in candidates if p.is_file()), None)


def canonical(key, aliases):
    seen = set()
    while key in aliases and aliases[key] != key:
        if key in seen:
            raise ValueError('Cycle in EVENT_ALIASES: ' + key)
        seen.add(key)
        key = aliases[key]
    return key


def add_structures(path, events):
    """Stream the evidence once; tolerate alias records arriving after contigs."""
    aliases, structures = {}, []
    with open_text(path) as stream:
        for number, line in enumerate(stream, 1):
            try:
                item = json.loads(line)
            except (ValueError, TypeError) as exc:
                raise ValueError(f'{path}:{number}: invalid JSON') from exc
            if item.get('record_kind') == 'EVENT_ALIASES':
                aliases.update(item['aliases'])
            elif item.get('record_kind') == 'CONTIG':
                contexts = {c['context_id']: c for c in item.get('contexts', [])}
                atoms = {a['atom_id']: a for a in item.get('atoms', [])}
                for event in item.get('events', []):
                    if event.get('svtype') not in ('CPX', 'DUP', 'DELINS'):
                        continue
                    members = list(event.get('components', [])) + [event.get('basis_atom')]
                    structures.append(dict(event=event, context=contexts.get(event.get('context_id'), {}),
                                           atoms=[atoms[k] for k in dict.fromkeys(members) if k in atoms]))
    for structure in structures:
        key = canonical(structure['event']['event_id'], aliases)
        if key in events:
            events[key]['structures'].append(structure)


def breakpoint_entries(event):
    """Component-coordinate endpoints, deliberately not a claim about wiring."""
    entries = []

    def add(kind, chrom, start, end, length):
        if any(c in chrom for c in '-:'):
            raise ValueError(f'VERIX single breakpoint syntax cannot encode contig {chrom!r}')
        if start < 0 or end < start:
            raise ValueError(f'Invalid component interval: {chrom}:{start}-{end}')
        cuts = [f'{chrom}:{start}']
        if end != start:
            cuts.append(f'{chrom}:{end}')
        entry = f'{kind}:{abs(length)}-' + '-'.join(cuts)
        if entry not in entries:
            entries.append(entry)

    info = event['main']['info']
    kind = event['main']['type']
    intervals = values(info.get('CPX_INTERVALS')) if kind == 'CPX' else []
    if kind == 'DUP':
        c = event['main']
        add('DUP', c['chrom'], c['pos'], c['end'], c['end']-c['pos'])
        if 'TARGET_POS' not in info or 'TARGET_CHR' not in info:
            raise ValueError(f"{event['id']}: native complex DUP lacks an explicit target")
    elif kind == 'DELINS':
        c = event['main']
        if int(info.get('INS_LEN', 0)) <= 0 or int(info.get('DEL_LEN', 0)) <= 0:
            raise ValueError(f"{event['id']}: DELINS requires INS_LEN and DEL_LEN")
        add('DEL', c['chrom'], c['pos'], c['end'], int(info['DEL_LEN']))
        add('INS', c['chrom'], c['pos'], c['pos'], int(info['INS_LEN']))
    elif not intervals:
        raise ValueError(f"{event['id']}: CPX_INTERVALS missing; cannot safely export the framework")
    for interval in intervals:
        match = re.fullmatch(r'([^_]+)_(.+):(\d+)-(\d+)', interval)
        if not match:
            raise ValueError(f'Invalid CPX_INTERVALS entry: {interval}')
        kind, chrom, first, last = match.groups()
        start, end = int(first) - 1, int(last)
        add(kind, chrom, start, end, end - start)
    for child in event['components']:
        ci = child['info']
        kind, chrom, start, end = child['type'], child['chrom'], child['pos'], child['end']
        if event['main']['type']=='DUP' and kind=='INS' and ci.get('INS_ORIGIN')=='DUPLICATION':
            # This is the same copied allele's INS view, not another novel INS.
            continue
        if kind == 'INS':
            add('INS', chrom, start, start, int(ci['SVLEN']))
        elif kind == 'DELINS':
            add('DEL', chrom, start, end, end - start)
            if 'INS_LEN' not in ci:
                raise ValueError(f"{child['id']}: DELINS component lacks INS_LEN")
            add('INS', chrom, start, start, int(ci['INS_LEN']))
        elif kind in ('DEL', 'INV', 'DUP'):
            add(kind, chrom, start, end, end - start)
        else:
            raise ValueError(f"{child['id']}: unsupported component type {kind}")
    if 'TARGET_POS' in info:
        pos = int(info['TARGET_POS'])
        add('TARGET', info.get('TARGET_CHR', event['main']['chrom']), pos, pos, 0)
    return entries


def load_events(vcf, components, evidence=None, selected=(), pass_only=False, scope='csv'):
    if scope not in ('csv', 'cpx'):
        raise ValueError('Unknown export scope: ' + scope)
    events = {}
    selected = set(selected)
    for fields, info, columns in records(vcf):
        category = csv_class(info, fields[4])
        if category is None or scope=='cpx' and info.get('SVTYPE')!='CPX' or selected and fields[2] not in selected:
            continue
        if pass_only and fields[6] != 'PASS':
            continue
        key = fields[2]
        if key == '.' or key in events:
            raise ValueError('Composite records require unique non-missing ID: ' + key)
        events[key] = dict(id=key, main=compact_record(fields, info, columns),
                           csv_class=category, components=[], structures=[], warnings=[], fields=fields)
    missing = selected - set(events)
    if missing:
        raise ValueError('Requested CSV absent, outside scope, or filtered: ' + ', '.join(sorted(missing)))
    if components is None and any(e['main']['type']=='CPX' for e in events.values()):
        raise ValueError('Component VCF not found; supply --components (expected <calls.vcf>.components.vcf)')
    if components:
        for fields, info, columns in records(components):
            parents = set(values(info.get('EVENT')) + values(info.get('PARENT_EVENT')))
            for parent in parents & events.keys():
                events[parent]['components'].append(compact_record(fields, info, columns))
                if len(parents) > 1:
                    events[parent]['warnings'].append('SHARED_COMPONENT')
    if evidence:
        LOG.info('Reading structural evidence: %s', evidence)
        add_structures(evidence, events)
    for event in events.values():
        if event['main']['type']=='CPX' and event['csv_class']=='COMPOSITE' and event['structures']:
            categories={composite_class(st['event'].get('info', {})) for st in event['structures']}
            if len(categories)==1:
                event['csv_class']=categories.pop()
        if not event['structures']:
            event['warnings'].append('NO_STRUCTURAL_EVIDENCE' if event['main']['type']=='CPX' else 'VCF_ONLY_LAYOUT')
        if not event['components'] and event['main']['type']=='CPX':
            event['warnings'].append('NO_COMPONENT_ROWS')
        event['warnings'] = sorted(set(event['warnings']))
        event['bkps'] = breakpoint_entries(event)
    return list(events.values())


def encode_parts(parts):
    return '|'.join(quote(str(p), safe='._-') for p in parts)


def write_vcf(path, source, events, scope='csv'):
    with open_text(source) as stream, open(path, 'w', encoding='utf-8') as out:
        event_declared = False
        for line in stream:
            if line.startswith('##INFO=<ID=EVENT,'):
                event_declared = True
            if any(line.startswith('##INFO=<ID=' + key + ',') for key in FIELDS):
                raise ValueError('Input already contains CPX export fields; use the original caller VCF')
            if line.startswith('#CHROM\t'):
                out.write('##asmsv_csv_export=' + SCHEMA + '\n')
                out.write('##asmsv_csv_scope=' + scope + '\n')
                out.write('##asmsv_csv_selection=' + POLICY + '\n')
                out.write('##asmsv_cpx_scope=Component-coordinate_breakpoints_not_full_ALT_or_junction_equivalence\n')
                if not event_declared:
                    out.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="Canonical parent event ID">\n')
                for key, (number, kind, description) in FIELDS.items():
                    out.write(f'##INFO=<ID={key},Number={number},Type={kind},Description="{description}">\n')
                out.write(line)
                break
            out.write(line)
        for event in events:
            info = parse_info(event['fields'][7])
            info.update(EVENT=event['id'], BKPS=','.join(event['bkps']), ASMSV_CSV_CLASS=event['csv_class'])
            children, genotypes = [], []
            for c in event['components']:
                ci = c['info']
                children.append(encode_parts([c['id'], c['type'], c['chrom'], c['pos'], c['end'],
                                              ci.get('SVLEN', '.'), ci.get('EDIT_CONTEXT', '.')]))
                for sample, data in c['samples'].items():
                    genotypes.append(encode_parts([c['id'], sample, data.get('GT', '.'),
                                                  data.get('PS', '.'), data.get('FT', '.')]))
            if children:
                info['ASMSV_COMPONENTS'] = ','.join(children)
            if genotypes:
                info['ASMSV_COMPONENT_GT'] = ','.join(genotypes)
            if event['warnings']:
                info['ASMSV_EXPORT_WARN'] = ','.join(event['warnings'])
            row = event['fields'].copy()
            row[7] = ';'.join(k if v is True else k + '=' + str(v) for k, v in info.items())
            out.write('\t'.join(row) + '\n')


def view_payload(events, provenance, scope='csv'):
    return dict(schema=SCHEMA, scope=scope, selection_policy=POLICY, inputs=provenance,
                events=[{k: v for k, v in e.items() if k != 'fields'} for e in events])


def write_jsonl(path, payload):
    with open(path, 'w', encoding='utf-8') as out:
        out.write(json.dumps(dict(record_kind='EXPORT_METADATA', schema=SCHEMA,
                                  scope=payload['scope'], selection_policy=POLICY,
                                  inputs=payload['inputs']), ensure_ascii=False) + '\n')
        for event in payload['events']:
            out.write(json.dumps(dict(record_kind='CSV_EVENT', **event), ensure_ascii=False) + '\n')


def write_html(path, payload):
    # Derived display models stay in HTML; canonical JSONL and VERIX stay intact.
    display = dict(payload, diagrams={e['id']: diagrams_for_event(e) for e in payload['events']})
    data = json.dumps(display, ensure_ascii=False, separators=(',', ':'))
    for char, escaped in [('&', '\\u0026'), ('<', '\\u003c'), ('>', '\\u003e'),
                          ('\u2028', '\\u2028'), ('\u2029', '\\u2029')]:
        data = data.replace(char, escaped)
    Path(path).write_text(HTML.replace('/*CPX_DATA*/', data), encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description='Export composite SV JSONL, VERIX single VCF and an offline HTML browser.')
    parser.add_argument('--vcf', required=True, type=Path, help='Original main calls.vcf[.gz]')
    parser.add_argument('--components', type=Path, help='Default: companion .components.vcf[.gz]')
    parser.add_argument('--events', type=Path, help='Default: companion .events.jsonl[.gz], if present')
    parser.add_argument('--no-events', action='store_true', help='Skip optional structural JSONL evidence')
    parser.add_argument('--out-prefix', required=True, type=Path, help='Example: calls.cpx')
    parser.add_argument('--scope', choices=('csv','cpx'), default='csv',
                        help='csv: CPX, DELINS and inverted/dispersed DUP; cpx: only SVTYPE=CPX')
    parser.add_argument('--event', action='append', default=[], help='Only this composite record ID; repeatable')
    parser.add_argument('--pass-only', action='store_true', help='Select PASS parent events only')
    parser.add_argument('--force', action='store_true', help='Replace existing export files; never inputs')
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
    temps = []
    try:
        if args.no_events and args.events:
            raise ValueError('--events and --no-events are mutually exclusive')
        components = args.components or companion(args.vcf, '.components.vcf')
        evidence = None if args.no_events else args.events or companion(args.vcf, '.events.jsonl')
        inputs = [p for p in (args.vcf, components, evidence) if p is not None]
        for path in inputs:
            if not path.is_file():
                raise ValueError('Input not found: ' + str(path))
        outputs = [Path(str(args.out_prefix) + suffix) for suffix in ('.verix.vcf', '.html', '.jsonl')]
        for path in outputs:
            if path.resolve() in {p.resolve() for p in inputs}:
                raise ValueError('Output cannot replace an input: ' + str(path))
            if path.exists() and not args.force:
                raise ValueError('Output exists; use --force to replace: ' + str(path))
        LOG.info('Reading composite events and linked components (scope=%s)', args.scope)
        events = load_events(args.vcf, components, evidence, args.event, args.pass_only, args.scope)
        provenance = {key: str(path.resolve()) if path else None for key, path in
                      [('main', args.vcf), ('components', components), ('events', evidence)]}
        provenance['main_sha256'] = hashlib.sha256(args.vcf.read_bytes()).hexdigest()
        for path in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
            os.close(fd)
            temps.append(Path(name))
        write_vcf(temps[0], args.vcf, events, args.scope)
        payload = view_payload(events, provenance, args.scope)
        write_html(temps[1], payload)
        write_jsonl(temps[2], payload)
        for temporary, path in zip(temps, outputs):
            os.replace(temporary, path)
        LOG.info('Exported %d composite events; %d have structural evidence', len(events),
                 sum(bool(e['structures']) for e in events))
        for path in outputs:
            LOG.info('Written: %s', path)
        return 0
    except (OSError, ValueError, KeyError) as exc:
        LOG.error('%s', exc)
        return 2
    finally:
        for path in temps:
            path.unlink(missing_ok=True)


HTML = r'''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AsmSV · CSV browser</title>
<style>
:root{font-family:system-ui,sans-serif;color:#182b42;background:#f4f6fa}*{box-sizing:border-box}body{margin:0}
header{background:#142c48;color:white;padding:18px 28px}h1{font-size:23px;margin:0 0 8px}p{line-height:1.65}
header p{margin:0;color:#d8e4f1;font-size:13px}.layout{display:grid;grid-template-columns:255px minmax(0,1fr);gap:20px;padding:22px}
aside{align-self:start;position:sticky;top:16px}input{width:100%;padding:12px;border:1px solid #b6c4d4;border-radius:8px}
#list{max-height:76vh;overflow:auto;margin-top:12px}button{font:inherit;cursor:pointer}#list button{display:block;width:100%;text-align:left;
padding:12px;background:white;border:1px solid #dbe2ec;border-radius:7px;margin:0 0 7px;overflow-wrap:anywhere;font-size:12px}
#list button.active{border:2px solid #227ba1;background:#effaff}article{min-width:0}.card{background:white;border:1px solid #dbe2ec;border-radius:10px;padding:20px;margin-bottom:16px}
h2{font-size:18px;margin:0 0 12px}h3{font-size:15px;margin:20px 0 10px}.mono,code{font-family:ui-monospace,monospace;overflow-wrap:anywhere}
.muted{font-size:13px;color:#52677d}.notice{background:#fff6df;border-left:4px solid #d79820;padding:10px 14px;line-height:1.7;font-size:13px}
.badge{display:inline-block;background:#e6f3f7;border-radius:5px;padding:5px 8px;margin:0 7px 7px 0;font-size:13px}
table{border-collapse:collapse;width:100%;font-size:12px}th,td{text-align:left;padding:10px 8px;border-bottom:1px solid #e4eaf1;overflow-wrap:anywhere}
.scroll{overflow-x:auto}svg{display:block;min-width:650px;width:100%;background:#fafcff;border-radius:8px;margin-bottom:12px;font-size:12px;fill:#233e58}svg text{font-family:system-ui,sans-serif}
pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:11px;max-height:500px;overflow:auto;background:#f5f7fa;padding:12px}summary{cursor:pointer}
.event-title{font-size:25px;font-weight:750;letter-spacing:-.4px;margin:9px 0 15px}.event-id{font-size:12px;color:#6b7b8d;font-weight:400}
.allele-panel{margin:20px -4px 8px;padding:16px 4px 8px;border-top:1px solid #e4eaf1}.allele-panel h3{margin:0 12px 6px;font-size:16px}.allele-panel>p{margin:0 12px 10px}
svg.alleleplot{min-width:1000px;background:white;margin:0;max-width:none}
.event-tree{border:1px solid #e0e7ee;background:#f8fafc;border-radius:9px;padding:14px 18px;margin:6px 10px 14px}
.event-tree ul{padding-left:22px;margin:8px 0}.event-tree li{padding:5px 0;font-size:14px;line-height:1.5}
.projection{display:flex;gap:12px;align-items:baseline;flex-wrap:wrap;border-top:1px solid #dee6ed;margin-top:10px;padding-top:12px;font-size:13px}
.projection strong{font-size:18px;color:#9c6417}.diagram-caption{font-size:12px;color:#64778b;margin:0 12px 12px}
.evidence>summary{font-weight:600;color:#52677d}.evidence>summary::marker{color:#8699ac}.evidence .card{box-shadow:none;border-color:#e7edf3}
@media(max-width:850px){.layout{grid-template-columns:1fr;padding:12px}aside{position:static}#list{max-height:200px}header{padding:18px}}
</style>
<header><h1>AsmSV · CSV 结构对照</h1><p>复合结构的 REF／ALT、来源与目标；保留原始 SVTYPE</p></header>
<div class="layout"><aside><input id="search" aria-label="搜索事件" placeholder="搜索 ID、染色体、类型或组件…"><p id="count" class="muted"></p><nav id="list"></nav></aside>
<article id="detail"></article></div>
<script id="data" type="application/json">/*CPX_DATA*/</script>
<script>
'use strict';
const data=JSON.parse(document.getElementById('data').textContent), all=data.events;
const $=id=>document.getElementById(id), NS='http://www.w3.org/2000/svg'; let current=all[0]?.id;
function node(tag,text,parent,cls){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;if(parent)parent.append(e);return e;}
function shape(svg,tag,attrs,text){const e=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;svg.append(e);return e;}
function canvas(parent,height){const box=node('div',undefined,parent,'scroll'),s=document.createElementNS(NS,'svg');s.setAttribute('viewBox',`0 0 940 ${height}`);s.setAttribute('role','img');box.append(s);return s;}
function table(parent,heads,rows){const box=node('div',undefined,parent,'scroll'),t=node('table',undefined,box),h=node('tr',undefined,node('thead',undefined,t));heads.forEach(x=>node('th',x,h));const b=node('tbody',undefined,t);rows.forEach(r=>{const tr=node('tr',undefined,b);r.forEach(x=>node('td',String(x??'.'),tr));});}
function formats(samples){return Object.entries(samples).map(([s,v])=>s+': '+Object.entries(v).map(([k,x])=>k+'='+x).join(' ')).join('；')||'无样本字段';}
const colors={INV:'#7659b4',DEL:'#c66565',DUP:'#398d9b',INS:'#cc933e'};
const contextNames={REFERENCE:'原位参考编辑',ADDED_DUPLICATE:'新增拷贝内编辑',DUPLICATION_JUNCTION:'复制连接处编辑'};
const classNames={COMPOSITE:'复合结构（未细分）',INV_WITH_EDITS:'倒位伴随编辑',DUP_WITH_EDITS:'复制伴随编辑',MULTISOURCE_INSERTION:'多来源插入',REPLACEMENT:'DEL＋INS 替换',INVERTED_DUPLICATION:'倒向复制',DISPERSED_DUPLICATION:'分散复制',INVERTED_DISPERSED_DUPLICATION:'倒向分散复制'};
function eventLabel(e){const i=e.main.info;if(e.main.type==='CPX')return 'CPX · '+(i.CPX_COMPONENTS||'').split(',').join(' ＋ ');if(e.main.type==='DELINS')return 'DELINS · DEL ＋ INS';return 'DUP · '+(i.DUP_TYPE||'')+(i.DUP_ORIENT?' · '+i.DUP_ORIENT:'');}
function renderList(){const q=$('search').value.toLowerCase();const filtered=all.filter(e=>JSON.stringify([e.id,e.csv_class,e.main.chrom,e.main.info,e.components.map(c=>c.id)]).toLowerCase().includes(q));$('count').textContent=`显示 ${filtered.length} / ${all.length} 个 ${data.scope==='cpx'?'CPX':'CSV'}`;$('list').replaceChildren();if(!filtered.some(e=>e.id===current)){current=filtered[0]?.id;renderDetail();}filtered.forEach(e=>{const b=node('button',e.main.chrom+':'+e.main.pos+' · '+eventLabel(e)+'\n'+e.id,$('list'),e.id===current?'active':'');b.onclick=()=>{current=e.id;renderList();renderDetail();};});}
function intervalPlot(parent,intervals){const rows=intervals.map(v=>{const m=/^([^_]+)_(.+):(\d+)-(\d+)$/.exec(v);return m?{kind:m[1],chr:m[2],a:Number(m[3]),b:Number(m[4])}:null;}).filter(Boolean);if(!rows.length)return;
 const s=canvas(parent,50+rows.length*58);shape(s,'text',{x:18,y:21},'CPX_INTERVALS · 1-based inclusive；此图不是 ALT 路径');
 rows.forEach((r,i)=>{const group=rows.filter(x=>x.chr===r.chr),a=Math.min(...group.map(x=>x.a)),b=Math.max(...group.map(x=>x.b)),scale=x=>300+(x-a)/Math.max(1,b-a)*590,y=51+i*58;
 shape(s,'text',{x:18,y:y+4},r.kind+' · '+r.chr+':'+r.a+'–'+r.b);shape(s,'line',{x1:300,y1:y,x2:890,y2:y,stroke:'#d8e1ec','stroke-width':2});
 shape(s,'rect',{x:scale(r.a),y:y-10,width:Math.max(3,scale(r.b)-scale(r.a)),height:20,rx:4,fill:colors[r.kind]||'#647e98'});});}
function alleleDiagram(parent,event,model,index){
 const panel=node('section',undefined,parent,'allele-panel');
 node('h3',`REF / ALT · ${model.haplotype||'结构示意'}`,panel);
 node('p',model.query_name||'未提供结构证据',panel,'muted');
 const wrap=node('div',undefined,panel,'scroll'),s=document.createElementNS(NS,'svg');
 s.setAttribute('viewBox','0 0 1160 415');s.setAttribute('class','alleleplot');s.setAttribute('role','img');
 s.setAttribute('aria-label','REF 与 ALT 的结构对照');s.setAttribute('font-family','system-ui,sans-serif');s.setAttribute('font-size','13');s.dataset.status=model.status;wrap.append(s);
 const text=(x,y,t,attrs={})=>shape(s,'text',{x,y,...attrs},t);
 const line=(x1,y1,x2,y2,attrs={})=>shape(s,'line',{x1,y1,x2,y2,stroke:'#a8b8c9','stroke-width':2,...attrs});
 const tip=(el,t)=>{const title=document.createElementNS(NS,'title');title.textContent=t;el.append(title);};
 const num=n=>Number(n).toLocaleString('en-US');
 function flank(label,x,y){shape(s,'rect',{x,y:y-19,width:52,height:38,rx:6,fill:'#edf1f5'});text(x+26,y+5,label,{'text-anchor':'middle','font-weight':700,fill:'#64758a'});}
 const ry=108,ay=284;
 text(14,ry-2,'REF',{'font-size':19,'font-weight':750,fill:'#294d7b'});text(14,ry+18,'参考',{'font-size':12,fill:'#7d8c9b'});
 text(14,ay-2,'ALT',{'font-size':19,'font-weight':750,fill:'#24745e'});text(14,ay+18,'变异',{'font-size':12,fill:'#7d8c9b'});
 for(const y of [ry,ay]){line(118,y,1105,y,{stroke:'#d3dce5','stroke-width':2.5});flank('L',105,y);flank('R',1083,y);}
 if(model.status==='supported'&&model.layout==='source_target'){
  // Draw two explicit loci. Their separation is not a physical-scale axis.
  s.replaceChildren();
  text(16,ry+5,'REF',{'font-size':19,'font-weight':750,fill:'#294d7b'});text(16,ay+5,'ALT',{'font-size':19,'font-weight':750,fill:'#24745e'});
  const[sc,ss,se]=model.source,[tc,tp]=model.target;
  text(330,35,'来源：'+sc+':'+num(ss+1)+'–'+num(se),{'text-anchor':'middle','font-weight':600});
  text(890,35,'目标：'+tc+':'+num(tp)+'（端口）',{'text-anchor':'middle','font-weight':600});
  const copyArrow=(x,y,reverse,added)=>{const w=230,head=17,d=reverse?`M${x+w},${y-23}H${x+head}L${x},${y}L${x+head},${y+23}H${x+w}Z`:`M${x},${y-23}H${x+w-head}L${x+w},${y}L${x+w-head},${y+23}H${x}Z`;
   shape(s,'path',{d,fill:added?'#e0f3e9':'#e3ecfa',stroke:added?'#287b58':'#386697','stroke-width':1.4,'data-role':added?'added':'retained','data-direction':reverse?'-':'+'});
   text(x+w/2,y+6,'S '+(reverse?'−':'+'),{'text-anchor':'middle','font-size':18,'font-weight':700,fill:added?'#287b58':'#386697'});};
  for(const y of[ry,ay]){line(130,y,535,y);flank('L',120,y);flank('R',520,y);line(690,y,1090,y);flank('L',680,y);flank('R',1080,y);}
  copyArrow(215,ry,false,false);copyArrow(215,ay,false,false);copyArrow(775,ay,model.copy_direction==='-',true);
  text(890,ry-20,'参考目标位置',{'text-anchor':'middle',fill:'#7d8c9b'});
  shape(s,'path',{d:`M445,${ry+24} C560,190 700,190 890,${ay-35}`,fill:'none',stroke:'#a9bfd8','stroke-width':1.6,'stroke-dasharray':'5 5'});
  text(330,ay+66,'来源保留',{'text-anchor':'middle',fill:'#386697'});
  text(890,ay-47,'DUP · '+(model.copy_direction==='-'?'反向新增拷贝 ←':'同向新增拷贝 →'),{'text-anchor':'middle',fill:'#287b58','font-weight':650});
  text(890,ay+66,num(se-ss)+' bp',{'text-anchor':'middle',fill:'#287b58'});
  text(110,407,'来源与目标分开示意；不按两地距离拉伸，也不将反向复制虚构为另一条独立 INV。',{'font-size':11,fill:'#8291a1'});
 }else if(model.status!=='supported'){
  shape(s,'rect',{x:250,y:ry-20,width:570,height:40,rx:5,fill:'#dce8fa'});text(535,ry+5,'已报告的参考框架',{'text-anchor':'middle'});
  shape(s,'rect',{x:250,y:ay-22,width:570,height:44,rx:5,fill:'#f8fafb',stroke:'#a4b2c0','stroke-dasharray':'6 5'});
  text(535,ay+5,'ALT 相对排布待核对',{'text-anchor':'middle',fill:'#718397'});
  text(200,365,'当前证据不足以可靠画出完整片段顺序；请展开下方证据查看已报告的区间与连接。',{fill:'#8b6a30'});
 }else{
  const weight=t=>t.kind==='missing'?.32:t.kind==='insertion'?1.5:1;
  const gw=g=>g.kind==='insertion'?1.6:Math.max(1.65,g.tokens.reduce((v,t)=>v+weight(t),0)+.25);
  const items=model.alt,total=items.reduce((v,g)=>v+gw(g),0),gap=22,available=890-gap*(items.length-1);
  let cursor=178;const placements=items.map(g=>{const p={g,x:cursor,w:available*gw(g)/total};cursor+=p.w+gap;return p;});
  const retained=placements.find(p=>p.g.role==='retained');
  const reference=model.ref[0],rp=retained?{g:reference,x:retained.x,w:retained.w}:{g:reference,x:178,w:890};
  const added=placements.find(p=>p.g.role==='added');
  if(added){shape(s,'path',{d:`M${rp.x+rp.w/2},${ry+24} C${rp.x+rp.w/2},190 ${added.x+added.w/2},185 ${added.x+added.w/2},${ay-46}`,
   fill:'none',stroke:'#b5c8df','stroke-width':1.8,'stroke-dasharray':'5 5'});text((rp.x+rp.w/2+added.x+added.w/2)/2,197,'复制来源',{'text-anchor':'middle',fill:'#7895b8','font-size':12});}
  function draw(p,y){const{g,x,w}=p;
   if(g.kind==='insertion'){
    shape(s,'rect',{x,y:y-22,width:w,height:44,rx:8,fill:'#fff0cf',stroke:'#d69b3c','stroke-width':1.5,'data-role':g.role||'insertion'});
    text(x+w/2,y+6,g.label,{'text-anchor':'middle','font-size':20,'font-weight':750,fill:'#9c661a'});
    text(x+w/2,y-44,'INS '+num(g.length)+' bp',{'text-anchor':'middle','font-weight':650,fill:'#aa6b0c'});
    text(x+w/2,y+56,g.title||'插入',{'text-anchor':'middle',fill:'#aa6b0c','font-size':12});
    return;
   }
   const role=g.role,bg=role==='added'?'#e0f3e9':role==='inverted'?'#ece5f7':'#e3ecfa';
   const ink=role==='added'?'#287b58':role==='inverted'?'#7855a6':'#386697';
   const groupNode=shape(s,'g',{'data-role':role});
   shape(groupNode,'rect',{x:x-7,y:y-30,width:w+14,height:60,rx:9,fill:bg,opacity:.38});
   const sum=g.tokens.reduce((v,t)=>v+weight(t),0),unit=(w-5*(g.tokens.length-1))/sum;
   let left=x;
   for(const t of g.tokens){const tw=unit*weight(t),cx=left+tw/2;
    if(t.kind==='missing'){
     line(cx-5,y-15,cx+5,y+15,{stroke:'#cb5c61','stroke-width':2.5});line(cx+5,y-15,cx-5,y+15,{stroke:'#cb5c61','stroke-width':2.5});
     text(cx,y-31,'Δ'+t.label,{'text-anchor':'middle',fill:'#b74d55','font-size':12});
     const mark=shape(groupNode,'rect',{x:left,y:y-23,width:tw,height:46,fill:'transparent','data-kind':'missing','data-label':t.label});
     tip(mark,'DEL '+t.length+' bp；'+(role==='added'?'新增拷贝缺少 ':'缺少 ')+t.label);
     line(cx,y+20,cx,y+46,{stroke:'#c96b70','stroke-width':1.2});
     text(cx,y+66,'DEL '+num(t.length)+' bp',{'text-anchor':'middle',fill:'#b74d55','font-size':12});
    }else if(t.kind==='insertion'){
     const el=shape(groupNode,'rect',{x:left,y:y-20,width:tw,height:40,rx:6,fill:'#fff0cf',stroke:'#d69b3c','data-kind':'insertion'});
     text(cx,y+5,t.label,{'text-anchor':'middle',fill:'#a66c16','font-weight':700});tip(el,'INS '+t.length+' bp；'+(contextNames[t.context]||'框架内插入'));
     text(cx,y+66,'INS '+num(t.length)+' bp',{'text-anchor':'middle',fill:'#a66c16','font-size':12});
    }else{
     const rev=t.direction==='-',head=Math.min(13,tw/5),d=rev?`M${left+tw},${y-20}H${left+head}L${left},${y}L${left+head},${y+20}H${left+tw}Z`:
      `M${left},${y-20}H${left+tw-head}L${left+tw},${y}L${left+tw-head},${y+20}H${left}Z`;
     const el=shape(groupNode,'path',{d,fill:bg,stroke:ink,'stroke-width':1.3,'data-kind':'segment','data-label':t.label,'data-direction':t.direction});
     text(cx,y+6,t.label+(t.direction==='-'?' −':' +'),{'text-anchor':'middle',fill:ink,'font-size':17,'font-weight':700});
     tip(el,`${t.label} ${t.direction} · ${t.chrom}:[${t.start},${t.end}) · ${t.length} bp`);
    }left+=tw+5;
   }
   let title=g.title;if(role==='added')title='DUP · 新增拷贝 '+(g.direction==='-'?'←':'→');
   if(role==='inverted')title='INV · 反向排列 ←';
   text(x+w/2,y-46,title,{'text-anchor':'middle',fill:ink,'font-weight':650});
   if(role==='reference'){
    const a=g.source; text(x+w/2,y+56,`${a[0]}:${num(a[1]+1)}–${num(a[2])}`,{'text-anchor':'middle',fill:'#6c8098','font-size':12});
   }else if(role==='added'){
    const n=g.tokens.filter(t=>t.kind==='missing');
    text(x+w/2,y+94,n.length?'新增拷贝缺少 '+n.map(t=>t.label).join('、'):'来源片段的额外拷贝',{'text-anchor':'middle',fill:n.length?'#b74d55':ink,'font-size':12});
   }
  }
  draw(rp,ry);placements.forEach(p=>draw(p,ay));
 }
 if(model.layout!=='source_target')text(110,407,'大 SV 结构示意 · 统一到参考正向 · 非比例，细小差异未展开；箭头为方向，× / Δ 为该次遍历缺失',{'font-size':11,fill:'#8291a1'});
 const tree=node('div',undefined,panel,'event-tree');node('strong','父事件：'+eventLabel(event),tree);
 for(const a of model.annotations||[])node('p',a,tree,'muted');
 const list=node('ul',undefined,tree);for(const c of model.summaries||[])node('li',`${c.type} ${num(c.length)} bp：${c.text}`,list);
 for(const p of model.projections||[]){const box=node('div',undefined,tree,'projection');node('span','indels 投影',box);node('strong','INS '+num(p.length)+' bp',box);node('span','整个新增等位结构；与组件分别查看',box);}
 if(model.status!=='supported')node('p','ALT 排布未补推；原始组件和连接仍保留在证据详情中。',tree,'muted');
}
function structurePlot(parent,st,index){const e=st.event,c=st.context,box=node('div',undefined,parent,'card');node('h3',`结构证据 ${index+1} · ${e.haplotype} · ${e.query_name}`,box);
 node('p','CONTEXT_ID='+e.context_id+'；target_strand='+(c.target_strand||'?')+'；解析状态='+(e.info.RESOLUTION||c.resolution||'?'),box,'muted');
 const segs=[...(e.segments||[])].sort((a,b)=>a.query[0]-b.query[0]);if(segs.length){node('p','来源在原始 query 中的访问顺序。箭头表示来源对参考的链方向；宽度为示意。未补齐小编辑或未比对序列。',box,'muted');const s=canvas(box,115);const w=Math.min(400,820/segs.length);segs.forEach((v,i)=>{const x=35+i*w,y=40,rev=v.source[3]==='-',d=rev?`M${x+w-12},${y-15}H${x+18}L${x},${y}L${x+18},${y+15}H${x+w-12}Z`:`M${x},${y-15}H${x+w-30}L${x+w-12},${y}L${x+w-30},${y+15}H${x}Z`;shape(s,'path',{d,fill:'#398d9b'});shape(s,'text',{x,y:77},'#'+(i+1)+' '+v.source[0]+':'+v.source[1]+'–'+v.source[2]+' '+v.source[3]);shape(s,'text',{x,y:98},'q:'+v.query[0]+'–'+v.query[1]);});
 table(box,['访问','来源 [start,end)','原始链','相对目标方向','query [start,end)'],segs.map((v,i)=>[i+1,v.source.slice(0,3).join(':'),v.source[3],c.target_strand?(v.source[3]===c.target_strand?'同向':'反向'):'未知',v.query.join('–')]));}
 node('p','target='+JSON.stringify(e.target)+'；reference_span='+JSON.stringify(e.reference_span)+'；query_span='+JSON.stringify(e.query_span)+'。这些 JSON 坐标为 0-based。',box,'mono muted');
 const edges=e.junctions||[];node('h3','已有连接配对 · '+edges.length+' 条',box);if(edges.length){const s=canvas(box,30+edges.length*67);edges.forEach((j,i)=>{const y=33+i*67,port=p=>p.chrom+':'+p.position+p.side;shape(s,'line',{x1:250,y1:y,x2:670,y2:y,stroke:'#7659b4','stroke-width':3});shape(s,'circle',{cx:250,cy:y,r:5,fill:'#7659b4'});shape(s,'circle',{cx:670,cy:y,r:5,fill:'#7659b4'});shape(s,'text',{x:20,y:y+4},port(j.left));shape(s,'text',{x:690,y:y+4},port(j.right));shape(s,'text',{x:250,y:y+23},'query support: '+JSON.stringify(j.query_support)+' · '+j.junction_id);});}else node('p','该证据中未保存连接；不根据区间自动补连。',box,'muted');
 if(c.issues?.length)node('p','context issues: '+c.issues.join(', '),box,'notice');}
function renderDetail(){const root=$('detail');root.replaceChildren();const e=all.find(x=>x.id===current);if(!e){node('p','没有匹配的 CPX 事件。',root,'card');return;}
 const main=e.main,card=node('section',undefined,root,'card');node('h2',e.id,card,'mono event-id');node('div',eventLabel(e),card,'event-title');[main.chrom+':'+main.pos+'–'+main.end,main.info.CPX_TYPE||'',classNames[e.csv_class]||'', 'FILTER='+main.filter].filter(Boolean).forEach(x=>node('span',x,card,'badge'));
 node('p','父事件原始分型：'+formats(main.samples),card,'mono muted');
 if(e.warnings.length)node('p','导出提示：'+e.warnings.join(', '),card,'notice');
 const models=data.diagrams[e.id]||[];if(models.length)models.forEach((m,i)=>alleleDiagram(card,e,m,i));else alleleDiagram(card,e,{status:'partial',summaries:[],projections:[]},0);
 const evidence=node('details',undefined,root,'card evidence');node('summary','展开组件、参考坐标与连接证据',evidence);
 node('h3','框架与编辑区间',evidence);intervalPlot(evidence,(main.info.CPX_INTERVALS||'').split(',').filter(Boolean));
 node('h3','组件归属 · '+e.components.length+' 行',evidence);table(evidence,['组件 ID','类型','VCF POS–END','长度','编辑位置','原始分型'],e.components.map(c=>[c.id,c.type,c.chrom+':'+c.pos+'–'+c.end,c.info.SVLEN||'.',(contextNames[c.info.EDIT_CONTEXT]||c.info.EDIT_CONTEXT||'?')+(c.shared_parents.length?'（共享组件）':''),formats(c.samples)]));
 node('p','INV／DUP 框架可能只在父行表达；组件不一定逐项对应 CPX_COMPONENTS。新增拷贝内 DEL 不代表原位来源缺失。',evidence,'muted');
 if(e.structures.length)e.structures.forEach((s,i)=>structurePlot(evidence,s,i));else node('p','没有关联到 events.jsonl 结构证据，仅展示 VCF 区间与组件。',evidence,'notice');
 const end=node('details',undefined,root,'card evidence');node('summary','展开 VERIX 断点与汇总原始字段',end);node('p','每个复合事件一行，保留原始 SVTYPE；single 模式读取 BKPS。断点为 0-based interbase；不包含连接配对、方向或完整 ALT 序列等价判定。',end,'muted');node('pre','BKPS='+e.bkps.join(','),end);
 const d=node('details',undefined,end);node('summary','查看汇总原始字段与输入来源',d);node('pre',JSON.stringify({inputs:data.inputs,event:e},null,2),d);}
$('search').addEventListener('input',renderList);renderList();renderDetail();
</script></html>'''


if __name__ == '__main__':
    raise SystemExit(main())
