import tempfile
from pathlib import Path
import unittest

import pysam

from asmsv_next.config import SVConfig
from asmsv_next.genotype_evidence import GenotypeEvidence, contig_evidence, summarize_genotypes
from asmsv_next.genotype import _merge_candidates
from asmsv_next.presets import resolve_preset
from asmsv_next.structures import SVCandidate
from asmsv_next.vcf import write_vcf
from test_core import block


def bundle(hap, name, spans, phase=(0, 2000)):
    return dict(haplotype=hap, genotyping=dict(reference=[['chr1', s, e, name, name+'-aln'] for s, e in spans],
                  phase=[['chr1', *phase, name]]))


def alt(kind='DEL', sources=('Hap1',), **info):
    c = SVCandidate(kind, 'chr1', 500, 100, 'q1' if sources == ('Hap1',) else 'q2', info, end=600)
    c.sources = set(sources)
    return c


class GenotypeEvidenceTests(unittest.TestCase):
    def setUp(self):
        SVConfig.CALLER_PRESET = resolve_preset('minimap2')
        SVConfig.MIN_MAPQ = 10; SVConfig.MIN_SV_SIZE = 30
        SVConfig.DETECTOR = 'structured'; SVConfig.GENOTYPING = 'evidence'; SVConfig.SAMPLE_NAME = 'SAMPLE'

    def tearDown(self):
        SVConfig.SAMPLE_NAME = 'SAMPLE'

    def test_alt_is_never_changed_to_reference_and_absence_needs_coverage(self):
        e = GenotypeEvidence([bundle('hap1','q1',[(0,2000)]), bundle('hap2','q2',[])])
        c=alt(); e.apply([c])
        self.assertEqual((c.genotype,c.format_data), ('1/.',{'FT':'NO_REFERENCE_SUPPORT'}))
        self.assertEqual(c.sources,{'Hap1'})

    def test_default_assembly_dosage_is_qualified_and_strict_mode_is_explicit(self):
        e=GenotypeEvidence([bundle('hap1','q1',[(0,2000)])])
        c=alt();e.apply([c],policy='assembly')
        self.assertEqual(c.genotype,'0/1')
        self.assertEqual(c.format_data,{'FT':'NO_REFERENCE_SUPPORT'})
        e.apply([c],policy='evidence');self.assertEqual(c.genotype,'1/.')
        stats=summarize_genotypes([c]);self.assertEqual(stats['unfiltered_fraction'],0)

    def test_reference_and_phase_require_the_same_physical_contig_pair(self):
        e = GenotypeEvidence([bundle('hap1','q1',[(0,2000)]), bundle('hap2','q2',[(0,2000)])])
        c=alt(); e.apply([c])
        self.assertEqual(c.genotype,'1|0'); self.assertEqual(c.format_data,{'PS':1})
        d=alt(sources=('Hap2',)); e.apply([d])
        self.assertEqual(d.genotype,'0|1'); self.assertEqual(d.format_data,{'PS':1})

    def test_overlapping_contigs_preclude_a_reference_assertion(self):
        e=GenotypeEvidence([bundle('hap1','q1',[(0,2000)]),bundle('hap2','q2',[(0,2000)]),
                            bundle('hap2','q3',[(580,900)])])
        c=alt(); e.apply([c])
        self.assertEqual(c.genotype,'1/.'); self.assertEqual(c.format_data['FT'],'AMBIGUOUS_COVERAGE')

    def test_gaps_are_not_bridged_and_deletion_span_must_be_callable(self):
        e=GenotypeEvidence([bundle('hap2','q2',[(0,550),(550,2000)])])
        c=alt(); e.apply([c]); self.assertEqual(c.genotype,'1/.')

    def test_duplication_uses_target_and_added_copy_component_is_conditional(self):
        e=GenotypeEvidence([bundle('hap2','q2',[(0,800)])])
        c=alt('DUP',TARGET_CHR='chr1',TARGET_POS=1500)
        e.apply([c]); self.assertEqual(c.genotype,'1/.')
        d=alt(ROLE='COMPONENT',ALLELE_CONTEXT='ADDED_DUPLICATE')
        e.apply([d]); self.assertEqual(d.format_data['FT'],'NONREFERENCE_COMPONENT')

    def test_two_observed_alts_remain_homozygous_without_a_reference_requirement(self):
        c=alt(sources=('Hap1','Hap2')); GenotypeEvidence([]).apply([c])
        self.assertEqual(c.genotype,'1/1'); self.assertEqual(c.format_data,{})

    def test_cluster_preserves_partial_filter_and_counts_physical_contigs(self):
        a,b,c=alt(),alt(),alt(sources=('Hap2',))
        for row in (a,b,c):row.filter='PARTIAL'
        result=_merge_candidates([a,b,c])
        self.assertEqual(result.filter,'PARTIAL')
        self.assertEqual((result.info['SUPPORT_H1'],result.info['SUPPORT_H2']),(1,1))

    def test_reference_chunks_reuse_index_and_stop_at_structural_indels_both_strands(self):
        for reverse in (False, True):
            b=block([(7,200),(1,50),(7,100),(2,100),(7,200)],is_reverse=reverse)
            rows=contig_evidence([b],{0},'q')['reference']
            self.assertEqual([(r[1],r[2]) for r in rows],[(1000,1200),(1200,1300),(1400,1600)])
            small=block([(7,200),(1,2),(7,100),(2,2),(7,200)],is_reverse=reverse)
            rows=contig_evidence([small],{0},'q')['reference']
            self.assertEqual([(r[1],r[2]) for r in rows],[(1000,1502)])

    def test_phase_sets_split_at_either_haplotype_contig_change(self):
        e=GenotypeEvidence([bundle('hap1','q1',[],phase=(0,2000)),
            bundle('hap2','q2',[(0,1000)],phase=(0,1000)),
            bundle('hap2','q3',[(1000,2000)],phase=(1000,2000))])
        a=alt(); b=alt(); b.pos,b.end=1400,1500
        e.apply([a,b])
        self.assertNotEqual(a.format_data['PS'],b.format_data['PS'])

    def test_format_roundtrip_sample_name_and_missing_rate(self):
        e=GenotypeEvidence([bundle('hap1','q1',[(0,2000)]),bundle('hap2','q2',[(0,1000)])])
        a=alt(); b=alt(); b.pos,b.end=1500,1600
        e.apply([a,b]); SVConfig.SAMPLE_NAME='HG002'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'calls.vcf'
            write_vcf([a,b],str(path),{'chr1':2000},None)
            with pysam.VariantFile(str(path)) as stream:
                rows=list(stream); self.assertEqual(list(stream.header.samples),['HG002'])
                self.assertEqual(rows[0].samples[0]['PS'],1)
                self.assertTrue(rows[0].samples[0].phased)
                self.assertEqual(rows[1].samples[0]['GT'],(1,None))
                self.assertFalse(rows[1].samples[0].phased)
        stats=summarize_genotypes([a,b]); self.assertEqual(stats['call_rate'],.5)


if __name__=='__main__':unittest.main()
