"""Explicit research-export scope for current AsmSV composite records.

This policy does not relabel the caller's SVTYPE or infer mutation history.
DELINS is a replacement stratum, not a claim of uniquely decomposed edits.
"""

POLICY = 'asmsv-csv-selection-v1'


def composite_class(info):
    """Normalize known frame semantics; retain a conservative fallback."""
    frame=info.get('FRAME')
    parts=str(info.get('CPX_COMPONENTS','')).split(',')
    sources=[v for v in str(info.get('CPX_INTERVALS','')).split(',')
             if v.startswith(('DUP_','INS_'))]
    if frame=='MULTISOURCE' or info.get('INS_ORIGIN')=='TEMPLATED' and len(sources)>=2:
        return 'MULTISOURCE_INSERTION'
    if frame=='INV' or 'INV' in parts and any(p in ('INS','DEL','DELINS') for p in parts):
        return 'INV_WITH_EDITS'
    if frame=='DUP' or len(sources)==1 and sources[0].startswith('DUP_') and any(p in ('INS','DEL','DELINS') for p in parts):
        return 'DUP_WITH_EDITS'
    return 'COMPOSITE'


def csv_class(info, alt=''):
    kind = info.get('SVTYPE')
    if kind == 'CPX':
        return composite_class(info)
    if kind == 'DELINS':
        return 'REPLACEMENT'
    if kind != 'DUP':
        return None
    subtype = info.get('DUP_TYPE', info.get('SUBTYPE', ''))
    if alt == '<DUP:TANDEM>' or subtype == 'DUP:TANDEM':
        return None
    if subtype == 'dDUP':
        return ('INVERTED_DISPERSED_DUPLICATION' if info.get('DUP_ORIENT')=='INVERTED'
                else 'DISPERSED_DUPLICATION')
    if subtype in ('piDUP_FR', 'piDUP_RF') or info.get('DUP_ORIENT') == 'INVERTED':
        return 'INVERTED_DUPLICATION'
    return None
