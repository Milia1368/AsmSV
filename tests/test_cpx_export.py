import gzip
import json
from pathlib import Path
import re
import tempfile
import unittest
from urllib.parse import unquote

import pysam

from asmsv_next.cpx_export import load_events, main, records, view_payload, write_html


HEADER = '''##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##ALT=<ID=CPX,Description="Composite allele">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=EVENT,Number=1,Type=String,Description="Parent">
##INFO=<ID=PARENT_EVENT,Number=.,Type=String,Description="Parents">
##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Intervals">
##INFO=<ID=CPX_COMPONENTS,Number=.,Type=String,Description="Types">
##INFO=<ID=EDIT_CONTEXT,Number=1,Type=String,Description="Edit context">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
'''
PARENT = 'chr1\t100\tp1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=400;SVLEN=-200;CPX_COMPONENTS=INV,DEL,DEL;CPX_INTERVALS=INV_chr1:201-300,DEL_chr1:101-200,DEL_chr1:301-400\tGT:PS\t0|1:9\n'
CHILD = 'chr1\t100\tc1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200;SVLEN=-100;EDIT_CONTEXT=REFERENCE;EVENT=p1\tGT:PS\t1|0:9\n'


class CpxExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.vcf = self.root / 'calls.vcf'
        self.components = self.root / 'calls.vcf.components.vcf'
        self.vcf.write_text(HEADER + PARENT)
        self.components.write_text(HEADER + CHILD)
        self.prefix = self.root / 'calls.cpx'

    def run_export(self, *extra):
        return main(['--vcf', str(self.vcf), '--out-prefix', str(self.prefix), *extra])

    def test_one_event_coordinates_and_original_genotypes(self):
        self.assertEqual(self.run_export(), 0)
        with pysam.VariantFile(str(self.prefix) + '.verix.vcf') as stream:
            rows = list(stream)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.id, 'p1')
        self.assertEqual(row.info['EVENT'], 'p1')
        self.assertEqual(row.samples['S']['GT'], (0, 1))
        self.assertTrue(row.samples['S'].phased)
        self.assertEqual(row.samples['S']['PS'], 9)
        self.assertEqual(row.info['BKPS'], ('INV:100-chr1:200-chr1:300',
                                          'DEL:100-chr1:100-chr1:200', 'DEL:100-chr1:300-chr1:400'))
        fields = row.info['ASMSV_COMPONENT_GT'][0].split('|')
        self.assertEqual([unquote(x) for x in fields], ['c1', 'S', '1|0', '9', '.'])
        payload = [json.loads(x) for x in Path(str(self.prefix) + '.jsonl').read_text().splitlines()]
        self.assertEqual(payload[0]['record_kind'], 'EXPORT_METADATA')
        self.assertEqual(payload[1]['components'][0]['samples']['S']['GT'], '1|0')

    def test_aliases_after_contig_and_two_haplotypes(self):
        path = Path(str(self.vcf) + '.events.jsonl')
        items = []
        for hap in ['hap1', 'hap2']:
            event = dict(event_id='internal_' + hap, svtype='CPX', context_id='ctx_' + hap,
                         haplotype=hap, query_name='q_' + hap, components=['a'])
            items.append(dict(record_kind='CONTIG', events=[event],
                              contexts=[dict(context_id='ctx_' + hap, issues=[])],
                              atoms=[dict(atom_id='a')]))
        items.append(dict(record_kind='EVENT_ALIASES', aliases={'internal_hap1': 'p1', 'internal_hap2': 'p1'}))
        path.write_text('\n'.join(json.dumps(x) for x in items))
        events = load_events(self.vcf, self.components, path)
        self.assertEqual(len(events), 1)
        self.assertEqual(len(events[0]['structures']), 2)
        self.assertEqual(events[0]['warnings'], [])

    def test_shared_parent_is_preserved_not_regrouped(self):
        second = PARENT.replace('\tp1\t', '\tp2\t')
        self.vcf.write_text(HEADER + PARENT + second)
        self.components.write_text(HEADER + CHILD.replace('EVENT=p1', 'PARENT_EVENT=p1,p2'))
        events = load_events(self.vcf, self.components)
        self.assertEqual(len(events), 2)
        for e in events:
            self.assertEqual(e['components'][0]['shared_parents'], ['p1', 'p2'])
            self.assertIn('SHARED_COMPONENT', e['warnings'])

    def test_missing_components_fails_without_partial_outputs(self):
        self.components.unlink()
        self.assertEqual(self.run_export(), 2)
        self.assertFalse(Path(str(self.prefix) + '.verix.vcf').exists())

    def test_existing_exports_require_force_and_inputs_are_protected(self):
        self.assertEqual(self.run_export(), 0)
        vcf_out = Path(str(self.prefix) + '.verix.vcf')
        before = vcf_out.read_bytes()
        self.assertEqual(self.run_export(), 2)
        self.assertEqual(vcf_out.read_bytes(), before)
        self.assertEqual(self.run_export('--force'), 0)
        self.vcf = vcf_out
        self.assertEqual(self.run_export('--force'), 2)
        self.assertEqual(vcf_out.read_bytes(), before)

    def test_missing_requested_event_does_not_export_everything(self):
        self.assertEqual(self.run_export('--event', 'not-found'), 2)
        self.assertFalse(Path(str(self.prefix) + '.html').exists())

    def test_empty_selection_and_gzipped_inputs(self):
        self.vcf.write_text(HEADER + PARENT.replace('\tPASS\t', '\t.\t'))
        gz = Path(str(self.vcf) + '.gz')
        with gzip.open(gz, 'wt') as stream:
            stream.write(self.vcf.read_text())
        self.vcf = gz
        self.assertEqual(self.run_export('--pass-only'), 0)
        self.assertEqual(len(list(records(str(self.prefix) + '.verix.vcf'))), 0)

    def test_browser_json_cannot_close_script_or_create_markup(self):
        events = load_events(self.vcf, self.components)
        malicious = '</script><img src=x onerror=alert(1)>'
        events[0]['main']['info']['RNAMES'] = malicious
        path = self.root / 'safe.html'
        write_html(path, view_payload(events, {}))
        text = path.read_text()
        self.assertNotIn(malicious, text)
        embedded = re.search(r'<script id="data" type="application/json">(.*?)</script>', text, re.S).group(1)
        self.assertEqual(json.loads(embedded)['events'][0]['main']['info']['RNAMES'], malicious)

    def test_missing_framework_refuses_incomplete_endpoint_export(self):
        self.vcf.write_text(HEADER + re.sub(r';CPX_INTERVALS=[^\t]+', '', PARENT))
        self.assertEqual(self.run_export(), 2)
        self.assertFalse(Path(str(self.prefix) + '.verix.vcf').exists())

    def test_added_copy_and_insertions_are_not_indel_projections(self):
        self.vcf.write_text(HEADER + PARENT.replace('INV_chr1:201-300,DEL_chr1:101-200,DEL_chr1:301-400', 'DUP_chr1:101-400'))
        self.components.write_text(HEADER + CHILD.replace('REFERENCE', 'ADDED_DUPLICATE') +
                                  'chr1\t400\tc2\tN\tNAAAA\t.\tPASS\tSVTYPE=INS;END=400;SVLEN=4;EVENT=p1\tGT\t1/.\n')
        # An unrelated projection file must never become another component.
        Path(str(self.vcf) + '.indels.vcf').write_text(HEADER + CHILD.replace('\tc1\t', '\tprojection\t'))
        events = load_events(self.vcf, self.components)
        self.assertEqual([c['id'] for c in events[0]['components']], ['c1', 'c2'])
        self.assertEqual(events[0]['components'][0]['info']['EDIT_CONTEXT'], 'ADDED_DUPLICATE')
        self.assertIn('INS:4-chr1:400', events[0]['bkps'])


if __name__ == '__main__':
    unittest.main()
