import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

import pysam
import test_structure as fixture

from asmsv_next.config import SVConfig
from asmsv_next.alignment import AlignedSegment
from asmsv_next.genotype import cluster_and_genotype
from asmsv_next.output_views import junction_candidates
from asmsv_next.region_output import write_regions
from asmsv_next.sequence_check import reverse_complement
from asmsv_next.structural_detection import detect_structures


class UnresolvedRegionTests(unittest.TestCase):
    def setUp(self):
        example = fixture.StructureTests()
        example.setUp()
        SVConfig.SEED_MAX_BASES=20000
        self.ref = example.ref

    def detect(self, query, chain, backbone, hap='hap1'):
        return detect_structures(chain, '+', backbone, 'q', query, fixture.Reference(self.ref), hap)

    def multisource(self, filler='', edited=False):
        ref = self.ref
        first = ref[8000:8200]+ref[8240:8500] if edited else ref[8000:8500]
        query = ref[:2000]+first+filler+reverse_complement(ref[9000:9500])+ref[10000:10500]+ref[2000:4000]
        first_body = [(7, 200), (2, 40), (7, 260)] if edited else [(7, 500)]
        middle = 2000+len(first)+len(filler)
        chain = [fixture.block(0, [(7, 2000)], 0, '+', len(query)),
                 fixture.block(2000, first_body, 8000, '+', len(query)),
                 fixture.block(middle, [(7, 500)], 9000, '-', len(query)),
                 fixture.block(middle+500, [(7, 500)], 10000, '+', len(query)),
                 fixture.block(middle+1000, [(7, 2000)], 2000, '+', len(query))]
        return query, chain

    def projection(self, result, hap='hap1'):
        return junction_candidates(dict(query='q', haplotype=hap, contexts=result.contexts, atoms=result.atoms), {})

    def test_three_sources_keep_four_connections_without_inventing_event(self):
        query, chain = self.multisource()
        result = self.detect(query, chain, {0, 4})
        self.assertEqual(result.calls, [])
        self.assertEqual(len(result.regions), 1)
        region = result.regions[0]
        self.assertEqual(region['query_span'], [2000, 3500])
        self.assertEqual(region['members']['block_indices'], [1, 2, 3])
        self.assertEqual(region['labels'], ['UNRESOLVED_REGION', 'CRR'])
        self.assertEqual(len(region['junctions']), 4)
        self.assertTrue(all(e['adjacency_status'] == 'OBSERVED' for e in region['junctions']))
        self.assertFalse(region['counts_as_variant'])
        self.assertIsNone(region['genotype'])
        records = self.projection(result)
        self.assertEqual(len(records), 8)
        ids = {c.id for c in records}
        self.assertTrue(all(c.mateid in ids and c.info['REGION_ID'] == region['region_id'] for c in records))
        self.assertTrue(all(c.filter == 'PARTIAL' for c in records))

    def test_open_branch_keeps_internal_connections_and_open_boundary(self):
        query, chain = self.multisource()
        query, chain = query[:-2000], chain[:-1]
        # Rebuild clipping for the shortened physical contig.
        chain = [fixture.block(b.qs, [(op, n) for op, n in b.cigar if op not in (4, 5)], b.ref_start, b.strand, len(query)) for b in chain]
        result = self.detect(query, chain, {0})
        self.assertEqual([c.type for c in result.calls], ['BND', 'BND'])
        region = result.regions[0]
        self.assertEqual(region['boundaries']['right']['state'], 'OPEN')
        self.assertEqual(region['query_span'], [2000, len(query)])
        self.assertEqual(len(region['junctions']), 3)
        self.assertEqual(len(self.projection(result)), 6)

    def test_double_gap_order_conflict_has_region_without_fabricated_breakends(self):
        ref = self.ref
        query = ref[:2000]+ref[2100:2200]+ref[2000:2100]+ref[2200:4300]
        chain = [fixture.block(0, [(7, 2000)], 0, '+', len(query)), fixture.block(2300, [(7, 2000)], 2300, '+', len(query))]
        result = self.detect(query, chain, {0, 1})
        region = result.regions[0]
        self.assertEqual(region['reasons'], ['SEED_ORDER_CONFLICT'])
        self.assertEqual(region['labels'], ['UNRESOLVED_REGION'])
        self.assertEqual(region['reason_categories'], ['SEQUENCE_AMBIGUITY'])
        self.assertEqual(region['junctions'], [])
        self.assertEqual(region['unretained_transitions'][0]['reason'], 'TWO_SIDED_RESIDUAL')
        self.assertTrue(region['sequence_checks'])
        self.assertEqual(result.calls, [])

    def test_assembly_gap_and_budget_are_diagnostics_without_crr(self):
        for ambiguous in (True, False):
            ref = self.ref
            if ambiguous:
                query = ref[:2000]+'N'*150+ref[2300:4300]
                start, expected, category = 2150, 'ASSEMBLY_GAP', 'ASSEMBLY_GAP'
            else:
                SVConfig.SEED_MAX_BASES = 100
                query = ref[:2000]+ref[2100:2200]+ref[2000:2100]+ref[2200:4300]
                start, expected, category = 2300, 'LARGE_RESIDUAL_BUDGET', 'COMPUTE_LIMIT'
            chain = [fixture.block(0, [(7, 2000)], 0, '+', len(query)), fixture.block(start, [(7, 2000)], 2300, '+', len(query))]
            region = self.detect(query, chain, {0, 1}).regions[0]
            self.assertIn(expected, region['reasons'])
            self.assertIn(category, region['reason_categories'])
            self.assertNotIn('CRR', region['labels'])
            self.assertEqual(region['junctions'], [])
        SVConfig.SEED_MAX_BASES = 20000

    def test_query_conflict_does_not_poison_other_certified_connections(self):
        query, chain = self.multisource()
        chain[2] = fixture.block(chain[2].qs-50, [(7, 550)], 9000, '-', len(query))
        result = self.detect(query, chain, {0, 4})
        region = result.regions[0]
        self.assertIn('QUERY_OWNERSHIP_CONFLICT', region['reasons'])
        self.assertTrue(any(b['reason'] == 'QUERY_OWNERSHIP_CONFLICT' for b in region['unretained_transitions']))
        self.assertTrue(region['junctions'])
        self.assertTrue(all(e['query_support'][0] <= e['query_support'][1] for e in region['junctions']))

    def test_original_atom_keeps_pass_and_ownership_with_region_link(self):
        query, chain = self.multisource(edited=True)
        result = self.detect(query, chain, {0, 4})
        self.assertEqual([c.type for c in result.calls], ['DEL'])
        atom = result.calls[0]
        region = result.regions[0]
        self.assertEqual(atom.filter, 'PASS')
        self.assertNotIn('PARENT_EVENT', atom.info)
        self.assertIn(atom.info['ATOM_ID'], region['members']['atom_ids'])
        self.assertEqual(atom.info['REGION_ID'], region['region_id'])
        self.assertEqual(result.cigar_scans, len(chain))

    def test_resolved_inversion_has_no_unresolved_region(self):
        ref = self.ref
        query = ref[:2000]+reverse_complement(ref[2000:3000])+ref[3000:5000]
        chain = [fixture.block(0, [(7, 2000)], 0, '+', len(query)), fixture.block(2000, [(7, 1000)], 2000, '-', len(query)),
                 fixture.block(3000, [(7, 2000)], 3000, '+', len(query))]
        result = self.detect(query, chain, {0, 2})
        self.assertEqual(result.regions, [])
        self.assertEqual([c.type for c in result.calls], ['INV'])
        self.assertNotIn('REGION_ID', result.calls[0].info)

    def test_junction_sequence_and_reverse_query_give_same_breakend_alleles(self):
        query, chain = self.multisource(filler='AGC')
        original = self.detect(query, chain, {0, 4})
        reversed_chain = [fixture.block(len(query)-b.qe, [(op, n) for op, n in b.cigar if op not in (4, 5)], b.ref_start,
                                       '-' if b.strand == '+' else '+', len(query)) for b in reversed(chain)]
        reverse = detect_structures(reversed_chain, '-', {0, 4}, 'q', reverse_complement(query), fixture.Reference(self.ref), 'hap1')
        def alleles(result):
            return sorted((c.chrom, c.pos, c.to_vcf_record(i)[4]) for i, c in enumerate(self.projection(result)))
        self.assertEqual(alleles(original), alleles(reverse))
        self.assertEqual(sum('BND_INSERTION' in c.info for c in self.projection(original)), 2)

    def test_region_ids_remain_distinct_across_haplotypes_while_bnds_merge(self):
        query, chain = self.multisource()
        first = self.detect(query, chain, {0, 4}, 'hap1')
        second = self.detect(query, chain, {0, 4}, 'hap2')
        ids = {first.regions[0]['region_id'], second.regions[0]['region_id']}
        self.assertEqual(len(ids), 2)
        merged = cluster_and_genotype(self.projection(first), self.projection(second, 'hap2'))
        self.assertEqual(len(merged), 8)
        self.assertTrue(all(set(c.info['REGION_ID'].split(',')) == ids for c in merged))

    def test_region_export_keeps_ports_exact_and_uses_one_base_bed_projection(self):
        query, chain = self.multisource()
        result = self.detect(query, chain, {0, 4})
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory/'hap.json').write_text(json.dumps(dict(regions=result.regions)))
            report = write_regions(directory, [], directory/'calls.vcf', {}, {'chr1': len(self.ref)})
            records = [json.loads(line) for line in Path(report['jsonl']).read_text().splitlines()]
            self.assertEqual(records[1]['vcf_records'], [])
            target = next(i for i in records[1]['reference_intervals'] if i['role'] == 'TARGET')
            self.assertEqual((target['start'], target['end']), (2000, 2000))
            rows = [line.split('\t') for line in Path(report['bed']).read_text().splitlines() if not line.startswith('#')]
            self.assertTrue(any(row[:3] == ['chr1', '2000', '2001'] and row[6] == 'TARGET' for row in rows))

    def test_cross_chromosome_members_keep_separate_reference_projections(self):
        query, chain = self.multisource()
        old = chain[2]
        chain[2] = AlignedSegment(dict(q_name='q', ref_name='chr2', ref_start=old.ref_start, ref_end=old.ref_end,
            cigar=old.cigar, is_reverse=True, mapq=60, query_length=len(query)))
        region = self.detect(query, chain, {0, 4}).regions[0]
        self.assertEqual({i['chrom'] for i in region['reference_intervals']}, {'chr1', 'chr2'})
        self.assertEqual(sum(e['left']['chrom'] != e['right']['chrom'] for e in region['junctions']), 2)
        self.assertTrue(all(i['end']-i['start'] <= 500 for i in region['reference_intervals']))

    def test_short_exact_support_and_unknown_junction_sequence_are_not_certified(self):
        query, chain = self.multisource(filler='N')
        region = self.detect(query, chain, {0, 4}).regions[0]
        self.assertTrue(any(t['reason'] == 'ASSEMBLY_GAP' for t in region['unretained_transitions']))
        self.assertTrue(all(e.get('inserted_sequence', '') != 'N' for e in region['junctions']))
        query, chain = self.multisource()
        old = chain[2]
        chain[2] = fixture.block(old.qs, [(7, 15)], old.ref_end-15, '-', len(query))
        region = self.detect(query, chain, {0, 4}).regions[0]
        self.assertEqual(sum(t['reason'] == 'INSUFFICIENT_EXACT_ANCHOR' for t in region['unretained_transitions']), 2)

    def test_repeated_reference_visits_have_distinct_occurrences_and_edges(self):
        ref = self.ref
        query = ref[:2000]+(ref[8000:8500]+reverse_complement(ref[9000:9500]))*2+ref[2000:4000]
        chain = [fixture.block(0, [(7, 2000)], 0, '+', len(query))]
        for qs, rs, strand in ((2000,8000,'+'),(2500,9000,'-'),(3000,8000,'+'),(3500,9000,'-')):
            chain.append(fixture.block(qs, [(7, 500)], rs, strand, len(query)))
        chain.append(fixture.block(4000, [(7, 2000)], 2000, '+', len(query)))
        region = self.detect(query, chain, {0, 5}).regions[0]
        self.assertEqual(len(region['junctions']), 5)
        self.assertEqual(len({e['junction_id'] for e in region['junctions']}), 5)
        self.assertEqual(len({o['occurrence_id'] for o in region['observed_path']}), 6)

    def test_unknown_vcf_region_reference_is_rejected_without_publishing(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            records = directory/'main.records.jsonl'
            records.write_text(json.dumps(dict(view='main', record_id='SV.1', info=dict(REGION_ID='missing')))+'\n')
            output = directory/'calls.vcf'
            with self.assertRaisesRegex(ValueError, 'absent unresolved regions'):
                write_regions(directory, [records], output, {}, {'chr1': len(self.ref)})
            self.assertFalse(Path(str(output)+'.regions.jsonl').exists())
            self.assertFalse(list(directory.glob('*.tmp')))

    def test_precollected_record_links_match_streamed_record_index(self):
        query, chain = self.multisource()
        region = self.detect(query, chain, {0, 4}).regions[0]
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory/'hap.json').write_text(json.dumps(dict(regions=[region])))
            records = directory/'main.records.jsonl'
            records.write_text(json.dumps(dict(view='main', record_id='known', info=dict(REGION_ID=region['region_id'])))+'\n')
            first = write_regions(directory, [records], directory/'first', {}, {'chr1': len(self.ref)})
            links={region['region_id']:[dict(view='main',record_id='known')]}
            second = write_regions(directory, [], directory/'second', {}, {'chr1': len(self.ref)}, record_links=links)
            self.assertEqual(Path(first['jsonl']).read_bytes(), Path(second['jsonl']).read_bytes())
            self.assertEqual(Path(first['bed']).read_bytes(), Path(second['bed']).read_bytes())
            self.assertEqual(len(links[region['region_id']]), 1)
            spool=directory/'regions.spool.jsonl'
            spool.write_text(json.dumps(region)+'\n')
            compact={region['region_id']:[('main','known')]}
            third=write_regions(directory,[],directory/'third',{}, {'chr1':len(self.ref)},record_links=compact,region_source=spool)
            self.assertEqual(Path(first['jsonl']).read_bytes(), Path(third['jsonl']).read_bytes())
            self.assertEqual(Path(first['bed']).read_bytes(), Path(third['bed']).read_bytes())

    def test_deletion_inside_anchor_or_at_incoming_port_is_not_certified(self):
        from asmsv_next.unresolved_regions import _exact_anchor
        block = fixture.block(0, [(7, 100), (2, 20), (7, 100)], 0, '+', 200)
        self.assertTrue(_exact_anchor(block, 100, True, 16))
        self.assertFalse(_exact_anchor(block, 100, False, 16))
        self.assertFalse(_exact_anchor(block, 108, True, 16))


class RegionCLITests(unittest.TestCase):
    def test_no_main_calls_still_export_regions_and_link_merged_bnd_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            rng = random.Random(73)
            reference = ''.join(rng.choice('ACGT') for _ in range(40000))
            query = reference[:6000]+reverse_complement(reference[22000:25000])+reference[28000:31000]+reverse_complement(reference[33000:36000])+reference[6000:12000]
            fasta = directory/'ref.fa'
            fasta.write_text('>chr1\n'+reference+'\n')
            pysam.faidx(str(fasta))
            header = dict(HD=dict(VN='1.6', SO='queryname'), SQ=[dict(SN='chr1', LN=len(reference))],
                          PG=[dict(ID='minimap2', PN='minimap2', VN='test', CL='minimap2 -ax asm5 --eqx')])
            hts = pysam.AlignmentHeader.from_dict(header)
            for hap in (1, 2):
                with pysam.AlignmentFile(str(directory/f'hap{hap}.bam'), 'wb', header=hts) as stream:
                    for n, (qs, size, rs, strand) in enumerate(((0,6000,0,'+'),(6000,3000,22000,'-'),(9000,3000,28000,'+'),(12000,3000,33000,'-'),(15000,6000,6000,'+'))):
                        clips = (qs, len(query)-qs-size) if strand == '+' else (len(query)-qs-size, qs)
                        read = pysam.AlignedSegment(hts)
                        read.query_name = 'complex'
                        read.reference_id = 0
                        read.reference_start = rs
                        read.flag = (2048 if n else 0) | (16 if strand == '-' else 0)
                        read.mapping_quality = 60
                        read.cigartuples = ([(4, clips[0])] if clips[0] else [])+[(7, size)]+([(4, clips[1])] if clips[1] else [])
                        read.query_sequence = reverse_complement(query) if strand == '-' else query
                        stream.write(read)
            output = directory/'calls.vcf'
            command = [sys.executable, '-B', '-m', 'asmsv_next', '--aligner', 'minimap2', '--bam1', str(directory/'hap1.bam'),
                       '--bam2', str(directory/'hap2.bam'), '--ref', str(fasta), '-t', '1', '-o', str(output), '--junction-vcf']
            root = Path(__file__).resolve().parents[1]
            result = subprocess.run(command, text=True, capture_output=True, timeout=30,
                                    env=dict(os.environ, PYTHONPATH=str(root/'src')))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(list(pysam.VariantFile(str(output))), [])
            records = [json.loads(line) for line in Path(str(output)+'.regions.jsonl').read_text().splitlines()]
            regions = records[1:]
            self.assertEqual(len(regions), 2)
            self.assertEqual({r['haplotype'] for r in regions}, {'hap1', 'hap2'})
            self.assertTrue(all(len(r['junctions']) == 4 and len(r['vcf_records']) == 8 for r in regions))
            bnds = list(pysam.VariantFile(str(output)+'.junctions.vcf'))
            self.assertEqual(len(bnds), 8)
            ids = {r['region_id'] for r in regions}
            self.assertTrue(all(set(row.info['REGION_ID']) == ids for row in bnds))
            manifest = json.loads(Path(str(output)+'.run.json').read_text())
            self.assertEqual(manifest['status'], 'COMPLETE')
            self.assertEqual(manifest['region_summary']['counts']['regions'], 2)
            self.assertIn('regions_bed', manifest['view_sha256'])
