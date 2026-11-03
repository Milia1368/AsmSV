"""VCF presentation policy. Never changes candidate detection or genotyping."""
import hashlib,json

BASE={'SVTYPE','END','SVLEN','SEQ','MATEID','RNAMES','RESOLUTION','REGION_ID'}
BIOLOGY=BASE|{'DUP_TYPE','DUP_ORIENT','CPX_TYPE','CPX_COMPONENTS','CPX_INTERVALS',
    'INS_ORIGIN','SOURCE_STATE','SOURCE_STATES','SOURCE_CHR','SOURCE_START','SOURCE_END','SOURCE_ORIENT',
    'TARGET_CHR','TARGET_POS','INS_LEN','DEL_LEN','EDIT_INS_LEN','EDIT_DEL_LEN',
    'JUNCTION_INS_LEN','EVENT','PARENT_EVENT','EDIT_CONTEXT'}
RENAMES={'SUBTYPE':'DUP_TYPE','PARENT_EVENT':'EVENT','ALLELE_CONTEXT':'EDIT_CONTEXT','SOURCE':'DETECTION_METHOD'}


def normalize_info(info,kind):
    result={RENAMES.get(key,key):value for key,value in info.items() if key!='PARENT_EVENT'}
    if info.get('PARENT_EVENT'):
        parents=str(info['PARENT_EVENT']).split(',')
        if len(parents)==1:result['EVENT']=parents[0]
        else:
            result.pop('EVENT',None);result['PARENT_EVENT']=info['PARENT_EVENT']
    if kind in ('CPX','DUP','INV'):
        for key in ('INS_LEN','DEL_LEN'):
            if key in result:result['EDIT_'+key]=result.pop(key)
    return result


def select_info(info,record,view='main',level='minimal'):
    kind=info['SVTYPE'];result=normalize_info(info,kind)
    if level=='full':return result
    if level!='minimal':raise ValueError('Unknown VCF annotation profile: '+level)
    role=info.get('ROLE');parent=info.get('PARENT_EVENT')
    if view=='indels':allowed=BASE|{'EVENT','PARENT_EVENT'}
    elif view=='junctions':allowed=BASE|{'EVENT','PARENT_EVENT'}
    elif view=='components':allowed=BASE|{'EVENT','PARENT_EVENT','EDIT_CONTEXT','INS_LEN','DEL_LEN',
        'INS_ORIGIN','SOURCE_STATE','SOURCE_CHR','SOURCE_START','SOURCE_END','SOURCE_ORIENT'}
    else:allowed=BIOLOGY
    result={key:value for key,value in result.items() if key in allowed}
    if result.get('EVENT')==record[2] and not parent:result.pop('EVENT',None)
    if result.get('CPX_TYPE')=='COMPLEX':result.pop('CPX_TYPE')
    if result.get('DUP_TYPE')=='DUP':result.pop('DUP_TYPE')
    if record[4]=='<DUP:TANDEM>':
        result.pop('DUP_TYPE',None);result.pop('DUP_ORIENT',None)
    elif result.get('DUP_TYPE') in ('piDUP_FR','piDUP_RF'):result.pop('DUP_ORIENT',None)
    if kind not in ('DUP','CPX'):
        if kind=='INS' and view in ('main','components') and info.get('DUP_ORIENT') and 'SOURCE_ORIENT' not in result:
            result['SOURCE_ORIENT']=info['DUP_ORIENT']
        result.pop('DUP_TYPE',None);result.pop('DUP_ORIENT',None)
    if kind!='INS' or view not in ('main','components'):
        for key in ('INS_ORIGIN','SOURCE_STATE','SOURCE_CHR','SOURCE_START','SOURCE_END','SOURCE_ORIENT'):result.pop(key,None)
    elif result.get('INS_ORIGIN') in (None,'UNKNOWN'):
        result.pop('INS_ORIGIN',None);result.pop('SOURCE_STATE',None)
    elif result.get('INS_ORIGIN')=='DUPLICATION':result.pop('SOURCE_STATE',None)
    # Multi-source insertion origin is meaningful on a CPX, independent of type.
    if kind=='CPX' and info.get('INS_ORIGIN')=='TEMPLATED':result['INS_ORIGIN']='TEMPLATED'
    if 'UNKNOWN' not in str(result.get('SOURCE_STATES','')).split(','):result.pop('SOURCE_STATES',None)
    if not (kind=='DUP' and info.get('DUP_MODE')!='TANDEM' or kind=='CPX' and 'DUP' in str(info.get('CPX_COMPONENTS','')).split(',')):
        result.pop('TARGET_CHR',None);result.pop('TARGET_POS',None)
    if kind!='DELINS':result.pop('INS_LEN',None);result.pop('DEL_LEN',None)
    for key in ('EDIT_INS_LEN','EDIT_DEL_LEN','JUNCTION_INS_LEN'):
        if str(result.get(key,'0'))=='0':result.pop(key,None)
    if result.get('RESOLUTION')=='RESOLVED' or record[6]=='PARTIAL' and result.get('RESOLUTION')=='PARTIAL':result.pop('RESOLUTION',None)
    if view!='components' and not (role=='COMPONENT' and view not in ('indels','junctions')):result.pop('EDIT_CONTEXT',None)
    return result


def compact_record(record,view='main',level='minimal'):
    info={}
    for token in record[7].split(';'):
        if token and token!='.':
            key,sep,value=token.partition('=');info[key]=value if sep else None
    selected=select_info(info,record,view,level)
    result=list(record)
    result[7]=';'.join(key if value is None else key+'='+str(value) for key,value in selected.items()) or '.'
    return tuple(result)


def header_definitions(definitions,level):
    result=[];seen=set()
    descriptions={
        'DUP_TYPE':'Named duplication subtype when not already specified by ALT',
        'EDIT_CONTEXT':'REFERENCE ADDED_DUPLICATE or DUPLICATION_JUNCTION; context of a component edit',
        'EVENT':'One associated canonical main-record ID; not a mutational-history claim',
        'DETECTION_METHOD':'Internal detection layer; diagnostic profile only',
        'CPX_TYPE':'Exact named complex template; full profile also exposes the unclassified fallback COMPLEX',
        'EDIT_INS_LEN':'Sum of inserted component edits; duplicated template gain is excluded',
        'EDIT_DEL_LEN':'Sum of deleted component edits; copy-local loss is not source-locus loss'}
    for key,number,kind,description in definitions:
        renamed=RENAMES.get(key,key)
        if level=='minimal' and renamed not in BIOLOGY:continue
        if renamed not in seen:
            result.append((renamed,1 if renamed=='EVENT' else number,kind,descriptions.get(renamed,description)));seen.add(renamed)
    for original in ('INS_LEN','DEL_LEN'):
        result.append(('EDIT_'+original,1,'Integer',descriptions['EDIT_'+original]))
    result.append(('PARENT_EVENT','.','String','Multiple distinct associated main-record IDs; EVENT is used for a unique parent'))
    return result


def record_evidence(candidate,record,view):
    """An exact ID index makes omitted diagnostics recoverable from JSONL."""
    info={key:value for key,value in candidate.info.items() if key!='SEQ'}
    sequence=candidate.info.get('SEQ','')
    return dict(record_kind='VCF_RECORD',view=view,record_id=record[2],svtype=candidate.type,
        chromosome=record[0],position=int(record[1]),end=candidate.end,length=candidate.sv_len,
        genotype=candidate.genotype,format_data=candidate.format_data.copy(),haplotype_support=sorted(candidate.sources),filter=record[6],info=info,sequence_length=len(sequence),
        sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        allele_sha256=hashlib.sha256(json.dumps([record[3],record[4]],separators=(',',':')).encode()).hexdigest())
