"""Read-only breakend projections of canonical structural evidence."""
from asmsv_next.structures import SVCandidate
from asmsv_next.structural_detection import junction,identifier


def junction_candidates(bundle,aliases):
    name=bundle['query'];haplotype=bundle['haplotype'];records=[];seen={}
    def emit(edge,context=None,parent=None,atom_id=None,partial=False):
        # Group complete adjacency pairs by exact geometry. Independent fuzzy
        # clustering of their two breakends can leave dangling MATEIDs.
        ports=sorted((edge['left'],edge['right']),key=lambda p:(p['chrom'],p['position'],p['side']))
        key=identifier('ADJ',ports,partial)
        if key in seen:
            for candidate in seen[key]:
                for field,value in [('CONTEXT_ID',context),('PARENT_EVENT',aliases.get(parent,parent)),('ATOM_ID',atom_id)]:
                    if value:candidate.info[field]=','.join(sorted((set(str(candidate.info.get(field,'')).split(','))|{value})-{''}))
            return
        seen[key]=[]
        for suffix,local,remote,mate in [('a',ports[0],ports[1],'b'),('b',ports[1],ports[0],'a')]:
            position=local['position']+(local['side']=='L')
            pos2=remote['position']+(remote['side']=='L')
            info=dict(SOURCE='STRUCTURE_VIEW',ROLE='JUNCTION',JUNCTION_ID=key,
                CHR2=remote['chrom'],POS2=max(1,pos2),
                STRANDS=('+' if local['side']=='R' else '-')+('+' if remote['side']=='L' else '-'),
                RESOLUTION='PARTIAL' if partial else 'RESOLVED')
            if context:info['CONTEXT_ID']=context
            if parent:info['PARENT_EVENT']=aliases.get(parent,parent)
            if atom_id:info['ATOM_ID']=atom_id
            if edge.get('query_support'):
                info['QUERY_START'],info['QUERY_END']=sorted(edge['query_support'])
            candidate=SVCandidate('BND',local['chrom'],max(1,position),0,name,info,
                filter_str='PARTIAL' if partial else 'PASS')
            candidate.id=key+'.'+suffix;candidate.mateid=key+'.'+mate;records.append(candidate)
            seen[key].append(candidate)
    for context in bundle['contexts']:
        for edge in context.get('junctions',[]):
            emit(edge,context['context_id'],context.get('event_id'),partial=context['resolution']!='RESOLVED')
    for atom in bundle['atoms']:
        if atom.get('parent_event') or atom['type'] not in ('INS','DEL','DELINS'):continue
        parts=[p for p in atom['parts'] if p[0] in (1,2)]
        if not parts:
            qs,qe=atom['query'];rs=min(p[0] for p in atom['reference_consumption']);re=max(p[1] for p in atom['reference_consumption'])
            parts=[[1 if atom['type']=='INS' else 2,qs,qe,rs,re]]
        for index,(_,qs,qe,rs,re) in enumerate(parts):
            left=dict(chrom=atom['chrom'],position=rs,side='R')
            right=dict(chrom=atom['chrom'],position=re,side='L')
            edge=junction(haplotype,name,atom['atom_id'],left,right,[qs,qe],str(index))
            emit(edge,atom_id=atom['atom_id'])
    return records
