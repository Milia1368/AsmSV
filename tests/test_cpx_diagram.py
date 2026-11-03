from copy import deepcopy
import unittest

from asmsv_next.cpx_diagram import diagram_for_structure


def atom(key, kind, start, end, query, context):
    return dict(atom_id=key, type=kind, chrom='chr1', pos=start, end=end,
                query=query, info=dict(ALLELE_CONTEXT=context))


def copy_structure():
    return dict(event=dict(event_id='e', haplotype='hap1', query_name='q', context_id='c',
                           reference_span=['chr1',100,600], query_span=[0,950],
                           info=dict(FRAME='DUP',DUP_MODE='ADJACENT',SOURCE_STATE='RETAINED',
                                     TARGET_CHR='chr1',TARGET_POS=600,DUP_OCCURRENCE=1,
                                     DUP_ORIENT='DIRECT',DUP_SEQUENCE_LEN=400),
                           segments=[dict(query=[0,500],source=['chr1',100,600,'+']),
                                     dict(query=[550,950],source=['chr1',100,600,'+'])]),
                context=dict(target_strand='+',source_intervals=[['chr1',100,600,'+']],
                             insertion_projections=[dict(length=450,target=['chr1',600])]),
                atoms=[atom('d','DEL',300,400,[750,750],'ADDED_DUPLICATE'),
                       atom('x','INS',600,600,[500,550],'DUPLICATION_JUNCTION')])


def visible(group):
    return [(t['label'],t['direction']) for t in group['tokens'] if t['kind']!='missing']


class CpxDiagramTests(unittest.TestCase):
    def test_deleted_piece_stays_in_reference_and_retained_copy(self):
        result=diagram_for_structure(copy_structure())
        self.assertEqual(result['status'],'supported')
        old, insertion, new = result['alt']
        self.assertEqual([t['label'] for t in result['ref'][0]['tokens']],['A','D','B'])
        self.assertEqual(visible(old),[('A','+'),('D','+'),('B','+')])
        self.assertEqual(visible(new),[('A','+'),('B','+')])
        self.assertEqual(new['tokens'][1]['kind'],'missing')
        self.assertEqual(insertion['role'],'junction')
        self.assertEqual(insertion['length'],50)
        self.assertEqual(result['projections'][0]['length'],450)

    def test_reverse_query_keeps_reference_relative_copy_order(self):
        s=copy_structure();s['context']['target_strand']='-'
        s['event']['info']['DUP_OCCURRENCE']=0
        s['event']['segments']=[dict(query=[0,400],source=['chr1',100,600,'-']),
                                dict(query=[450,950],source=['chr1',100,600,'-'])]
        s['atoms'][0]['query']=[200,200];s['atoms'][1]['query']=[400,450]
        result=diagram_for_structure(s)
        self.assertEqual(result['status'],'supported')
        self.assertEqual([g['role'] for g in result['alt']],['retained','junction','added'])
        self.assertEqual(visible(result['alt'][-1]),[('A','+'),('B','+')])

    def test_target_before_source_places_added_copy_first(self):
        s=copy_structure();s['event']['info'].update(TARGET_POS=100,DUP_OCCURRENCE=0)
        s['event']['segments']=[dict(query=[0,400],source=['chr1',100,600,'+']),
                                dict(query=[450,950],source=['chr1',100,600,'+'])]
        s['atoms'][0]['query']=[200,200];s['atoms'][1]['query']=[400,450]
        result=diagram_for_structure(s)
        self.assertEqual(result['status'],'supported')
        self.assertEqual([g['role'] for g in result['alt']],['added','junction','retained'])

    def test_inverted_added_copy_reverses_piece_order_and_direction(self):
        s=copy_structure();s['event']['info']['DUP_ORIENT']='INVERTED'
        s['event']['segments'][1]['source'][3]='-'
        result=diagram_for_structure(s)
        self.assertEqual(result['status'],'supported')
        self.assertEqual(visible(result['alt'][-1]),[('B','-'),('A','-')])

    def test_unknown_source_or_conflicting_target_does_not_invent_alt(self):
        for update in [dict(SOURCE_STATE='UNKNOWN'),dict(TARGET_POS=900),dict(DUP_OCCURRENCE=0)]:
            s=copy_structure();s['event']['info'].update(update)
            result=diagram_for_structure(s)
            self.assertEqual(result['status'],'partial')
            self.assertIn('reason',result)

    def test_inversion_with_flanking_deletions(self):
        s=copy_structure();s['event']['info']=dict(FRAME='INV')
        s['context']['source_intervals']=[['chr1',200,400,'-']]
        s['event']['segments']=[dict(query=[1000,1200],source=['chr1',200,400,'-'])]
        s['event']['query_span']=[1000,1200]
        s['atoms']=[atom('d1','DEL',100,200,[1000,1000],'REFERENCE'),
                    atom('d2','DEL',400,600,[1200,1200],'REFERENCE')]
        result=diagram_for_structure(s)
        self.assertEqual(result['status'],'supported')
        self.assertEqual([t['label'] for t in result['ref'][0]['tokens']],['D1','A','D2'])
        self.assertEqual(visible(result['alt'][0]),[('A','-')])

    def test_internal_inversion_deletion_is_removed_from_reversed_path(self):
        s=copy_structure();s['event']['info']=dict(FRAME='INV')
        s['context']['source_intervals']=[['chr1',100,600,'-']]
        s['event']['segments']=[dict(query=[0,400],source=['chr1',100,600,'-'])]
        s['event']['query_span']=[0,400]
        s['atoms']=[atom('d','DEL',300,400,[200,200],'REFERENCE')]
        result=diagram_for_structure(s)
        self.assertEqual(result['status'],'supported')
        self.assertEqual(visible(result['alt'][0]),[('B','-'),('A','-')])

    def test_schematic_does_not_modify_original_evidence(self):
        s=copy_structure();original=deepcopy(s)
        diagram_for_structure(s)
        self.assertEqual(s,original)

    def test_two_insertions_at_one_cut_follow_query_order(self):
        s=copy_structure()
        s['atoms'].extend([atom('late','INS',250,250,[710,730],'ADDED_DUPLICATE'),
                           atom('early','INS',250,250,[700,710],'ADDED_DUPLICATE')])
        result=diagram_for_structure(s)
        self.assertEqual(result['status'],'supported')
        inserted=[t for t in result['alt'][-1]['tokens'] if t['kind']=='insertion']
        self.assertEqual([t['edit_id'] for t in inserted],['early','late'])


if __name__=='__main__':
    unittest.main()
