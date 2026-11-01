from dataclasses import asdict
import importlib.util
import math
from pathlib import Path
import pickle
import random
import sys
import unittest

from asmsv_next.aq import RawAlignment,EvidenceRoot
from asmsv_next.presets import resolve_preset
from asmsv_next.alignment import AlignedSegment
from asmsv_next.config import SVConfig
from asmsv_next.structural_detection import map_port,mapped_parts
from asmsv_next.alignment_tags import SourcePreset

SPEC=importlib.util.spec_from_file_location('aq_before_packing',Path(__file__).parent/'fixtures/aq_before_packing.py')
ORACLE=importlib.util.module_from_spec(SPEC);sys.modules[SPEC.name]=ORACLE;SPEC.loader.exec_module(ORACLE)


def exact(value):
    if isinstance(value,float):return value.hex()
    if isinstance(value,dict):return {key:exact(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return tuple(exact(item) for item in value)
    return value


class PackedIndexTests(unittest.TestCase):
    def test_random_sliced_scores_and_views_are_bitwise_equal_to_frozen_index(self):
        rng=random.Random(11012026)
        for mode in ('asm5','asm10','asm20'):
            preset=resolve_preset('minimap2',mode).aq
            for case in range(40):
                body=[(7,100)]
                for _ in range(25):body.extend([(rng.choice((1,2,8)),rng.randrange(1,400)),(7,rng.randrange(1,300))])
                body.append((7,100));cigar=[(5,37),*body,(5,53)]
                if case%3==0:cigar=[(0 if op in (7,8) else op,n) for op,n in cigar]
                nm=sum(n for op,n in body if op in (1,2,8))
                data=dict(q_name='q',ref_name='chr1',ref_start=1000,cigar=cigar,mapq=60,
                          is_reverse=bool(case%2),tags={'NM':nm})
                new=EvidenceRoot(RawAlignment.from_dict(data),preset)
                old=ORACLE.EvidenceRoot(ORACLE.RawAlignment.from_dict(data),preset)
                self.assertEqual(asdict(new.raw),asdict(old.raw))
                self.assertEqual(new.issues,old.issues)
                self.assertEqual([r._asdict() for r in new.index.runs],[asdict(r) for r in old.index.runs])
                for _ in range(30):
                    start,end=sorted(rng.sample(range(new.raw.qs,new.raw.qe+1),2))
                    self.assertEqual(exact(asdict(new.index.score(start,end))),exact(asdict(old.index.score(start,end))))
                    self.assertEqual(new.mapped_bounds(start,end),old.mapped_bounds(start,end))
                    if new.mapped_bounds(start,end) is not None:
                        self.assertEqual(new.cigar_view(start,end),old.cigar_view(start,end))

    def test_slotted_evidence_pickle_round_trip(self):
        raw=RawAlignment.from_dict(dict(q_name='q',ref_name='chr1',ref_start=0,cigar=[(7,500),(2,60),(7,500)],mapq=60))
        root=EvidenceRoot(raw,resolve_preset('minimap2').aq)
        copy=pickle.loads(pickle.dumps(root))
        self.assertEqual(asdict(copy.raw),asdict(raw))
        self.assertEqual(copy.index.score(raw.qs,raw.qe),root.index.score(raw.qs,raw.qe))
        self.assertFalse(hasattr(root.index.runs[0],'__dict__'))

    def test_shared_source_metadata_preserves_dictionary_input_contract(self):
        source=SourcePreset(aligner='minimap2',mode='asm5',overrides=('-A','2','-B','6'))
        data=dict(q_name='q',ref_name='chr1',ref_start=0,cigar=[(7,1000)],mapq=60)
        direct=RawAlignment.from_dict(dict(data,source=source))
        dictionary=RawAlignment.from_dict(dict(data,source=asdict(source)))
        self.assertIs(direct.source,source)
        self.assertEqual(asdict(direct),asdict(dictionary))

    def test_indexed_ports_and_local_parts_match_full_scan_at_indel_boundaries(self):
        SVConfig.CALLER_PRESET=resolve_preset('minimap2');SVConfig.MIN_MAPQ=10
        rng=random.Random(44219)
        def scan(block,query=None,reference=None,exact=False):
            q0,q1=query or (block.qs,block.qe)
            q0,q1=max(q0,block.qs),min(q1,block.qe)
            for run in block.score_root.index.runs:
                if run.op not in ((7,) if exact else (0,7,8)):continue
                s,e=max(q0,run.qs),min(q1,run.qe)
                if s>=e:continue
                if block.strand=='+':a,b=run.rs+s-run.qs,run.rs+e-run.qs
                else:a,b=run.re-(e-run.qs),run.re-(s-run.qs)
                if reference is not None:
                    lo,hi=max(a,reference[0]),min(b,reference[1])
                    if lo>=hi:continue
                    if block.strand=='+':s,e=s+lo-a,e-(b-hi)
                    else:s,e=s+b-hi,e-(lo-a)
                    a,b=lo,hi
                yield s,e,a,b
        for reverse in (False,True):
            for _ in range(30):
                cigar=[(5,13),(7,100)]
                for _ in range(15):cigar.extend([(rng.choice((1,2,8)),rng.randrange(1,100)),(7,rng.randrange(1,100))])
                cigar.extend([(7,100),(5,17)])
                block=AlignedSegment(dict(q_name='q',ref_name='chr1',ref_start=1000,
                    ref_end=1000+sum(n for op,n in cigar if op in (0,2,3,7,8)),cigar=cigar,mapq=60,is_reverse=reverse))
                points={block.qs-1,block.qe+1,block.qs,block.qe}
                points.update(x for r in block.score_root.index.runs for x in (r.qs,r.qe,(r.qs+r.qe)//2))
                for point in points:
                    expected=next((a+point-s if block.strand=='+' else b-(point-s)
                                   for s,e,a,b in scan(block) if s<=point<=e),None)
                    self.assertEqual(map_port(block,point),expected)
                for _ in range(15):
                    query=sorted(rng.sample(range(0,block.qe+30),2))
                    for exact_only in (False,True):
                        reference=(1100,1400) if rng.randrange(2) else None
                        self.assertEqual(list(mapped_parts(block,query,reference,exact_only)),list(scan(block,query,reference,exact_only)))


if __name__=='__main__':unittest.main()
