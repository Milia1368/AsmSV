from dataclasses import replace
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import os
import sys
import tempfile
from types import SimpleNamespace
import unittest
import pysam
from asmsv_next.alignment import AlignedSegment
from asmsv_next.config import SVConfig
from fixtures.legacy_detection import resolve_chain_svs
from asmsv_next.dp import solve_path,transition_cost,edge_terms,path_score
from asmsv_next.input_io import bam_groups,inspect_headers
from asmsv_next.presets import AQPreset,DPPreset,resolve_preset
from asmsv_next.structures import SVCandidate
from asmsv_next.constants import fixed_policy

ROOT=Path(__file__).resolve().parents[1]


def block(cigar,**changes):
    data=dict(q_name='q',ref_name='chr1',ref_start=1000,cigar=cigar,mapq=60,tags={})
    data.update(changes)
    data.setdefault('ref_end',data['ref_start']+sum(n for op,n in cigar if op in (0,2,3,7,8)))
    return AlignedSegment(data)


def node(qs,qe,rs,score,chrom='chr1',strand='+'):
    return SimpleNamespace(qs=qs,qe=qe,ref_name=chrom,ref_start=rs,ref_end=rs+qe-qs,
        as_score=score,strand=strand,cigar=((7,qe-qs),),is_valid=True,is_backbone=False)


class AQTests(unittest.TestCase):
    def setUp(self):SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_MAPQ=10

    def test_shared_scoring_keeps_source_identity_without_mode_coupling(self):
        presets=[resolve_preset(tool,mode) for tool,modes in (
            ('minimap2',('asm5','asm10','asm20')),('winnowmap2',('asm5','asm10','asm20')),
            ('vacmap',('asm',))) for mode in modes]
        self.assertTrue(all((p.aq,p.dp)==(AQPreset(),DPPreset()) for p in presets))
        self.assertEqual(len({p.scoring_fingerprint for p in presets}),1)
        self.assertEqual(len({p.fingerprint for p in presets}),len(presets))
        self.assertEqual([(p.aligner,p.mode) for p in presets][-1],('vacmap','asm'))
        self.assertEqual(presets[3].aligner,'winnowmap')
        with self.assertRaises(ValueError):resolve_preset('vacmap','asm20')

    def test_continuous_gap_boundary_and_monotonicity(self):
        p=SVConfig.CALLER_PRESET.aq
        self.assertEqual((p.gap_cost(49),p.gap_cost(50)),(130,131))
        self.assertLess(abs(p.gap_cost(50+1e-7)-p.gap_cost(50-1e-7)),1e-5)
        values=[p.gap_cost(i) for i in range(1,10001)]
        self.assertTrue(all(a<=b for a,b in zip(values,values[1:])))

    def test_mapq_once_after_repeated_view(self):
        a=block([(5,100),(7,1000),(5,200)],mapq=10)
        x=a.split_at_query_coords(300,800);y=x.split_at_query_coords(300,800)
        self.assertEqual((x.as_score,y.as_score,x.aq_sequence,x.aq),(400,400,1,.8))

    def test_negative_AS_structural_TP_can_be_rescued(self):
        a=block([(7,600),(2,3331),(7,600)],raw_as=-900,tags={'AS':-900,'NM':3331})
        self.assertTrue(a.is_valid);self.assertEqual(a.raw_as,-900)

    def test_raw_AS_is_diagnostic_and_never_substitutes_for_missing_differences(self):
        scored=[]
        for native in (-100000,0,100000,None,float('nan'),'malformed'):
            a=block([(7,500),(8,3),(2,213),(7,500)],raw_as=native)
            scored.append((a.aq,a.as_score,a.is_valid,a.summary))
        self.assertTrue(all(x==scored[0] for x in scored))
        self.assertFalse(block([(0,1000)],raw_as=100000).is_valid)

    def test_MAPQ_admission_and_weight_are_separate_fixed_policies(self):
        for quality,weight,valid in ((9,.8,False),(10,.8,True),(19,.8,True),(20,1.,True),(60,1.,True),(255,None,False)):
            a=block([(7,1000)],mapq=quality)
            self.assertEqual((a.mapq_factor,a.is_valid),(weight,valid))
            if weight is not None:self.assertEqual((a.aq,a.as_score),(weight,1000*weight))

    def test_low_positive_AQ_and_short_structural_support_are_not_hard_gated(self):
        for cigar in ([(7,100),(8,9),(7,100)],[(7,100),(2,213),(7,238)]):
            self.assertTrue(block(cigar).is_valid)

    def test_refine_keeps_its_original_stable_view(self):
        a=block([(7,12),(1,65),(7,1000)])
        self.assertEqual((a.qs,a.qe,a.ref_start),(77,1077,1012))
        self.assertFalse(a.parse_intra_alignment_sv(None,30,3))

    def test_only_nine_scoring_fields_remain_and_removed_options_fail(self):
        from dataclasses import fields
        self.assertEqual(len(fields(AQPreset))+len(fields(DPPreset)),9)
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/'p.json'
            for data in ({'schema_version':4},{'admission':{'policy':'known'}},
                         {'aq':{'min_aq':.2}},{'aq':{'rescue_min_score':0}},
                         {'dp':{'short_rule':'none'}},{'dp':{'topology_unit':'aq_gap'}}):
                f.write_text(json.dumps(data))
                with self.assertRaises(ValueError):resolve_preset('minimap2',preset_file=f)

    def test_common_override_retains_source_identity_and_fixes_numeric_evaluation(self):
        numeric=[]
        for tool in ('minimap2','winnowmap','vacmap'):
            p=resolve_preset(tool,preset_file=ROOT/'presets/common-assembly-hq-v5.json')
            self.assertEqual(p.aligner,tool);numeric.append((p.aq,p.dp))
        self.assertEqual(numeric,[numeric[0]]*3)


    def test_negative_reward_and_low_aq_reach_dp_with_a_small_positive_floor(self):
        a=block([(7,100),(8,20),(7,100)],raw_as=-180)
        self.assertTrue(a.is_valid)
        self.assertEqual(a.as_score,1)
        self.assertLess(a.aq,0)
        self.assertTrue(a.reward_rescued)
        self.assertIn('REWARD_FLOORED',a.quality_flags)
        self.assertEqual(solve_path([a],220).blocks,(a,))
        # No gain from a weak, fully overlapping alternative mapping.
        strong=block([(7,219),(5,1)])
        self.assertEqual(solve_path([strong,a],220).blocks,(strong,))

    def test_rescue_does_not_invent_difference_evidence_or_ignore_bad_tags(self):
        self.assertFalse(block([(0,1000)]).is_valid)
        self.assertFalse(block([(7,1000)],tags={'NM':-1}).is_valid)
        self.assertFalse(block([(7,1000)],mapq=255).is_valid)

    def test_unknown_M_and_NM_only_interval_precision(self):
        self.assertFalse(block([(0,1000)]).is_valid)
        a=block([(0,1000)],tags={'NM':10})
        self.assertTrue(a.is_valid);self.assertFalse(a.summary.positional)
        self.assertFalse(a.split_at_query_coords(100,900).is_valid)

    def test_reverse_view_uses_real_coordinates(self):
        a=block([(5,100),(7,100),(2,60),(7,100),(5,200)],is_reverse=True)
        x=a.split_at_query_coords(300,400)
        self.assertEqual((x.ref_start,x.ref_end,x.qs,x.qe,x.as_score),(1000,1100,300,400,100))

    def test_custom_preset_roundtrip_and_invalid_values(self):
        p=resolve_preset('vacmap')
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/'preset.json';f.write_text(json.dumps(p.to_dict()))
            self.assertEqual(resolve_preset('vacmap',preset_file=f).fingerprint,p.fingerprint)
            f.write_text(json.dumps({'aq':{'structural_scale':12}}))
            partial=resolve_preset('vacmap',preset_file=f)
            self.assertEqual((partial.aq.mismatch,partial.aq.structural_scale),(19,12))
            f.write_text(json.dumps({'aq':{'mismatch':19,'gap_open':[39,81],'gap_extend':[3,1]}}))
            self.assertEqual(resolve_preset('vacmap',preset_file=f).scoring_fingerprint,p.scoring_fingerprint)
            f.write_text(json.dumps({'aq':{'mismatch':2,'gap_open':[2,12],'gap_extend':[1,.5]}}))
            custom=resolve_preset('vacmap',preset_file=f)
            self.assertEqual((custom.aq.mismatch,custom.aq.gap_open,custom.aq.gap_extend),(2,(2,12),(1,.5)))
            self.assertNotEqual(custom.scoring_fingerprint,p.scoring_fingerprint)
        for kwargs in ({'mismatch':float('nan')},{'structural_gap_min':0},{'gap_open':()}):
            with self.assertRaises(ValueError):AQPreset(**kwargs)
        with self.assertRaises(ValueError):DPPreset(query_gap=-1)


class DPTests(unittest.TestCase):
    def setUp(self):SVConfig.CALLER_PRESET=resolve_preset('minimap2')

    def test_independent_exhaustive_objective_200_graphs(self):
        rng=random.Random(8143);p=SVConfig.CALLER_PRESET.dp
        def score(path,length):
            if not path:return -.5*length
            value=sum(x.as_score for x in path)-.5*(path[0].qs+length-path[-1].qe)
            for a,b in zip(path,path[1:]):
                gap=b.qs-a.qe
                value-=.25*max(0,gap)+max(0,-gap)*max(a.as_score/(a.qe-a.qs),b.as_score/(b.qe-b.qs))+.5*max(0,-gap-10000)
                if a.ref_name!=b.ref_name:value-=500+100*min(1.,max(0.,2-min(a.as_score,b.as_score)/500))
                if a.strand!=b.strand:
                    value-=200
                    if a.ref_name==b.ref_name:
                        exit_a=a.ref_end if a.strand=='+' else a.ref_start
                        entry_b=b.ref_start if b.strand=='+' else b.ref_end
                        value-=min(500,.1*max(0,abs(entry_b-exit_a)-max(a.ref_end-a.ref_start,b.ref_end-b.ref_start)))
                elif a.ref_name==b.ref_name:
                    gr=b.ref_start-a.ref_end if b.strand=='+' else a.ref_start-b.ref_end
                    value-=min(500,.1*abs(gap-gr))
            return value
        for _ in range(200):
            nodes=sorted([node(qs:=rng.randrange(0,150),qs+rng.randrange(10,900),rng.randrange(3000),
                rng.randrange(1,1500),rng.choice(('chr1','chr2')),rng.choice(('+','-'))) for i in range(7)],key=lambda n:(n.qs,n.qe))
            best=-.5*1500
            for mask in range(1<<len(nodes)):
                path=[n for i,n in enumerate(nodes) if mask>>i&1]
                if any(not (a.qs<b.qs and a.qe<b.qe) for a,b in zip(path,path[1:])):continue
                best=max(best,score(path,1500))
            result=solve_path(nodes,1500)
            self.assertAlmostEqual(result.score,best,places=7)
            self.assertAlmostEqual(score(result.blocks,1500),best,places=7)

    def test_empty_path_and_overlap_topology_cost(self):
        self.assertEqual(solve_path([],100).score,-50)
        a=node(0,20000,0,18000);b=node(5000,30000,100000,22000,'chr2')
        self.assertAlmostEqual(transition_cost(a,b),15000*.9+2500+500)


    def test_information_smoothing_preserves_weak_node_prior_and_uses_boundary_budgets(self):
        a=node(0,3000,0,3000);p=SVConfig.CALLER_PRESET.dp
        def cost(reward,span=1500):
            return edge_terms(a,node(3000,3000+span,0,reward,'chr2'),p)['chromosome_change']
        self.assertEqual([cost(r) for r in (1,420,500,750,1000,2000)],[600,600,600,550,500,500])
        self.assertEqual(cost(420,840),cost(420,1500))
        for boundary in (500,1000):
            self.assertLess(abs(cost(boundary-1e-6)-cost(boundary+1e-6)),1e-6)
        # Explicitly tie the information budgets to the configured prior.
        q=replace(p,cross_chromosome=250)
        self.assertEqual(edge_terms(a,node(3000,4000,0,375,'chr2'),q)['chromosome_change'],275)


    def test_topology_prior_and_each_boundary_are_charged(self):
        a=node(0,1000,1000,1000)
        col=node(1000,2000,2000,1000)
        inv=node(1000,2000,2000,1000,strand='-')
        tra=node(1000,2000,2000,1000,chrom='chr2')
        self.assertLess(transition_cost(a,col),transition_cost(a,inv))
        self.assertLess(transition_cost(a,inv),transition_cost(a,tra))
        both=node(1000,2000,2000,1000,chrom='chr2',strand='-')
        self.assertEqual(transition_cost(a,both),700)
        c=node(2000,3000,3000,1000)
        self.assertEqual(edge_terms(inv,c)['orientation_change'],200)
        self.assertEqual(sum(edge_terms(x,y)['orientation_change'] for x,y in zip([a,inv,c],[inv,c])),400)

    def test_local_overlap_uses_actual_evidence_instead_of_whole_block_density(self):
        a=block([(7,800),(8,20),(7,180),(5,500)])
        b=block([(5,500),(7,1000)],ref_start=1500)
        self.assertAlmostEqual(edge_terms(a,b)['shared_query'],500)
        # The later block has low quality in the shared interval, high quality outside.
        c=block([(5,500),(7,100),(8,20),(7,880)],ref_start=1500)
        local_a=a.score_root.index.score(500,1000).node_reward
        local_c=c.score_root.index.score(500,1000).node_reward
        self.assertAlmostEqual(edge_terms(a,c)['shared_query'],max(0,local_a,local_c))
        self.assertNotAlmostEqual(edge_terms(a,c)['shared_query'],500*c.density)

    def test_reverse_complement_preserves_objective_and_unique_path(self):
        nodes=[node(0,1800,1000,1700),node(1600,2200,8000,500,'chr2','-'),
               node(2100,4000,3000,1750)]
        def reverse(b):
            return SimpleNamespace(**{**vars(b),'qs':4000-b.qe,'qe':4000-b.qs,
                'strand':'-' if b.strand=='+' else '+'})
        first=solve_path(nodes,4000);second=solve_path([reverse(b) for b in nodes],4000)
        self.assertAlmostEqual(first.score,second.score)
        self.assertEqual([(b.ref_name,b.ref_start) for b in first.blocks],
                         [(b.ref_name,b.ref_start) for b in reversed(second.blocks)])
        self.assertAlmostEqual(path_score(nodes,4000),path_score([reverse(b) for b in reversed(nodes)],4000))

    def test_reverse_complement_with_real_eqx_and_deletion_overlap(self):
        a=block([(7,800),(2,60),(7,200),(5,1000)])
        b=block([(5,500),(7,100),(8,20),(7,980),(5,400)],ref_start=1500)
        def reverse(x):
            return AlignedSegment(dict(q_name=x.q_name,ref_name=x.ref_name,ref_start=x.ref_start,
                ref_end=x.ref_end,cigar=x.cigar,is_reverse=True,mapq=x.mapq))
        ar,br=reverse(a),reverse(b)
        self.assertAlmostEqual(transition_cost(a,b),transition_cost(br,ar))
        self.assertAlmostEqual(solve_path([a,b],2000).score,solve_path([ar,br],2000).score)


    def test_reference_cost_has_no_ten_kilobase_drop(self):
        a=node(0,1000,1000,1000)
        values=[transition_cost(a,node(1000,2000,2000+x,1000)) for x in (4999,5000,9999,10000,10001)]
        self.assertEqual(values,sorted(values));self.assertEqual(values[-1],500)

    def test_record_order_does_not_change_selected_path(self):
        nodes=[node(0,1000,1000,900),node(900,2000,1900,1000),node(950,1950,80000,800,'chr2')]
        chosen=[tuple((b.qs,b.qe,b.ref_name) for b in solve_path(order,2000).blocks) for order in itertools.permutations(nodes)]
        self.assertEqual(len(set(chosen)),1)


class LegacyTests(unittest.TestCase):
    def setUp(self):SVConfig.CALLER_PRESET=resolve_preset('minimap2')

    def test_conservative_replacement_gate(self):
        a=block([(7,1000),(5,4000)])
        for gq,gr,expected in ((1000,1000,[]),(1500,1000,['CSV']),(0,900,['DEL'])):
            qs=1000+gq;b=block([(5,qs),(7,1000),(5,4000-qs)],ref_start=2000+gr)
            a.is_backbone=b.is_backbone=True
            calls=resolve_chain_svs([a,b],'+',{0,1},'q','A'*5000)
            self.assertEqual([c.type for c in calls],expected)
            if expected==['CSV']:
                self.assertEqual((calls[0].info['INS_LEN'],calls[0].info['DEL_LEN']),(1500,1000))

    def test_reference_fetch_returns_sequence_not_boolean(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'r.fa';path.write_text('>chr1\n'+'ACGT'*100+'\n');pysam.faidx(str(path))
            with pysam.FastaFile(str(path)) as ref:
                c=SVCandidate('DEL','chr1',20,50,'q')
                self.assertEqual(c._safe_fetch_ref(ref,'chr1',20,70),'ACGT'*12+'AC')
                row=c.to_vcf_record(1,reference=ref)
                self.assertTrue(set(row[3])<=set('ACGT'))


def fixture(directory,aligner,hardclip=False):
    ref=directory/'ref.fa';ref.write_text('>chr1\n'+'A'*10000+'\n');pysam.faidx(str(ref))
    aligned='A'*1000+'C'*60+'A'*2000;whole='G'*1000+aligned+'T'*1000
    query=directory/'query.fa';query.write_text('>toy\n'+whole+'\n');pysam.faidx(str(query))
    mode='asm' if aligner=='vacmap' else 'asm5'
    command='vacmap -mode asm --H --fakecigar' if aligner=='vacmap' else f'{aligner} -ax asm5 --eqx'
    header={'HD':{'VN':'1.6','SO':'queryname'},'SQ':[{'SN':'chr1','LN':10000}],
        'PG':[{'ID':aligner,'PN':aligner,'VN':'test','CL':command}]}
    hts=pysam.AlignmentHeader.from_dict(header)
    a=pysam.AlignedSegment(hts);a.query_name='toy';a.reference_id=0;a.reference_start=1000
    a.mapping_quality=60;a.flag=0
    clip=5 if hardclip else 4
    a.cigartuples=[(clip,1000),(7,1000),(1,60),(7,1000),(2,200),(7,1000),(clip,1000)]
    a.query_sequence=aligned if hardclip else whole;a.set_tag('NM',260)
    if aligner!='vacmap':a.set_tag('AS',2578)
    for hap in (1,2):
        with pysam.AlignmentFile(str(directory/f'hap{hap}.bam'),'wb',header=hts) as out:
            if hap==1:out.write(a)
    return ref,query


class InputAndCLITests(unittest.TestCase):

    def test_different_upstream_score_overrides_are_provenance_without_a_preset_file(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);ref,_=fixture(d,'minimap2')
            for hap in (1,2):
                path=d/f'hap{hap}.bam'
                with pysam.AlignmentFile(str(path),'rb') as stream:
                    header=stream.header.to_dict();reads=list(stream)
                flags='-A2 -B6 -O4,24 -E2,1' if hap==1 else '-A3 -B15 -O6,36 -E3,1.5'
                header['PG'][0]['CL']='minimap2 -ax asm10 --eqx '+flags
                with pysam.AlignmentFile(str(path),'wb',header=header) as out:
                    for r in reads:out.write(r)
            target=d/'result.vcf'
            result=subprocess.run(['bash',str(ROOT/'run.sh'),'--aligner','minimap2',
                '--bam1',str(d/'hap1.bam'),'--bam2',str(d/'hap2.bam'),'--ref',str(ref),'-t','1','-o',str(target)],
                text=True,capture_output=True,timeout=45,env=dict(os.environ,ASMSV_CALL_PY=sys.executable))
            self.assertEqual(result.returncode,0,result.stderr)
            manifest=json.loads(Path(str(target)+'.run.json').read_text())
            self.assertEqual((manifest['preset']['mode'],manifest['mode_origin']),('asm10','header'))
            self.assertEqual(manifest['preset']['aq']['mismatch'],19)
            self.assertEqual(manifest['preset']['aq']['gap_open'],[39,81])
            self.assertEqual(manifest['preset']['aq']['gap_extend'],[3,1])
            self.assertIn('-B6',manifest['source_headers'][0]['overrides'])
            self.assertIn('-B15',manifest['source_headers'][1]['overrides'])
            self.assertEqual(manifest['scoring_policy'],'fixed-common')
            self.assertFalse(manifest['upstream_scoring_inherited'])
            self.assertEqual(manifest['scoring_sha256'],resolve_preset('vacmap').scoring_fingerprint)

    def test_source_aware_end_to_end_three_tools(self):
        for tool in ('minimap2','winnowmap','vacmap'):
            with self.subTest(tool=tool),tempfile.TemporaryDirectory() as t:
                d=Path(t);ref,query=fixture(d,tool,tool=='vacmap');out=d/'result.vcf';trace=d/'trace'
                cmd=['bash',str(ROOT/'run.sh'),'--aligner',tool,'--bam1',str(d/'hap1.bam'),'--bam2',str(d/'hap2.bam'),
                    '--ref',str(ref),'-t','2','-o',str(out),'--trace-dir',str(trace)]
                if tool=='vacmap':cmd+=['--query-fasta1',str(query)]
                run=subprocess.run(cmd,text=True,capture_output=True,timeout=45,env=dict(os.environ,ASMSV_CALL_PY=sys.executable))
                self.assertEqual(run.returncode,0,run.stderr)
                with pysam.VariantFile(str(out)) as stream:rows=list(stream)
                self.assertEqual(sorted((r.info['SVTYPE'],abs(r.info['SVLEN'])) for r in rows),[('DEL',200),('INS',60)])
                manifest=json.loads(Path(str(out)+'.run.json').read_text())
                self.assertEqual(manifest['status'],'COMPLETE')
                self.assertEqual(manifest['fixed_policy'],fixed_policy())
                self.assertEqual(manifest['preset']['schema_version'],5)
                detail=json.loads(next(trace.glob('*.json')).read_text())
                self.assertEqual(detail['blocks'][0]['source']['mode'],'asm' if tool=='vacmap' else 'asm5')
                self.assertEqual(manifest['preset']['aq']['mismatch'],19)

    def test_hard_clip_requires_original_sequence_and_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);ref,query=fixture(d,'vacmap',True);SVConfig.CALLER_PRESET=resolve_preset('vacmap')
            with self.assertRaisesRegex(ValueError,'no complete SEQ'):list(bam_groups(d/'hap1.bam','hap1'))
            query.write_text('>toy\n'+'C'*5060+'\n');Path(str(query)+'.fai').unlink();pysam.faidx(str(query))
            with self.assertRaisesRegex(ValueError,'sequence endpoint sample differs'):list(bam_groups(d/'hap1.bam','hap1',str(query)))

    def test_declared_aligner_conflict_is_an_error(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);ref,query=fixture(d,'minimap2')
            with pysam.FastaFile(str(ref)) as fa:
                with self.assertRaisesRegex(ValueError,'conflicts'):inspect_headers((d/'hap1.bam',),fa,'vacmap')

    def test_text_CG_is_rejected_instead_of_dropped(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);ref,query=fixture(d,'vacmap',True)
            with pysam.AlignmentFile(str(d/'hap1.bam'),'rb') as bam:header=bam.header
            a=pysam.AlignedSegment(header);a.query_name='bad';a.flag=0;a.reference_id=0
            a.reference_start=100;a.mapping_quality=60;a.query_sequence='A'*100;a.set_tag('CG','100=')
            path=d/'text_cg.bam'
            with pysam.AlignmentFile(str(path),'wb',header=header) as out:out.write(a)
            SVConfig.CALLER_PRESET=resolve_preset('vacmap')
            with self.assertRaisesRegex(ValueError,'text CG:Z'):list(bam_groups(path,'hap1',str(query)))

    def test_noncontiguous_query_groups_are_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);ref,query=fixture(d,'minimap2')
            with pysam.AlignmentFile(str(d/'hap1.bam'),'rb') as bam:a=next(bam);header=bam.header
            path=d/'bad_order.bam'
            with pysam.AlignmentFile(str(path),'wb',header=header) as out:
                for name in ('toy','other','toy'):
                    a.query_name=name;out.write(a)
            SVConfig.CALLER_PRESET=resolve_preset('minimap2')
            with self.assertRaisesRegex(ValueError,'not contiguous'):list(bam_groups(path,'hap1'))


if __name__=='__main__':unittest.main()
