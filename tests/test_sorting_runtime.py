import json
import hashlib
import os
from pathlib import Path
import tempfile
import subprocess
import unittest
from unittest.mock import patch

import pysam

from asmsv_next.cli import identity, run_stage, sorted_bam


class SortingRuntimeTests(unittest.TestCase):
    def make_bam(self, directory):
        path = directory / 'input.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'chr1', 'LN': 1000}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as stream:
            for index, name in enumerate(('q2', 'q1', 'q2')):
                row = pysam.AlignedSegment(stream.header)
                row.query_name = name
                row.query_sequence = 'ACGT' * 10
                row.reference_id = 0
                row.reference_start = index * 50
                row.mapping_quality = 60
                row.cigarstring = '40M'
                stream.write(row)
        return path

    def test_sort_uses_separate_scratch_and_cache_is_reused(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            work = directory / 'cache'; work.mkdir()
            scratch = directory / 'scratch'; scratch.mkdir()
            source = self.make_bam(directory)
            original = subprocess.run
            with patch('asmsv_next.cli.subprocess.run', wraps=original) as sorting:
                result = sorted_bam(source, None, work, 1, scratch)
                argv = sorting.call_args.args[0]
                self.assertTrue(Path(argv[argv.index('-T') + 1]).is_relative_to(scratch))
            with pysam.AlignmentFile(result) as stream:
                self.assertEqual([row.query_name for row in stream], ['q1', 'q2', 'q2'])
            self.assertEqual(list(scratch.iterdir()), [])
            with patch('asmsv_next.cli.subprocess.run', side_effect=AssertionError('cache miss')):
                self.assertEqual(sorted_bam(source, None, work, 1, scratch), result)

    def test_default_caches_live_beside_each_input_and_preserve_sources(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            for hap in ('hap1', 'hap2'):
                parent = directory / hap; parent.mkdir()
                source = self.make_bam(parent)
                before = hashlib.sha256(source.read_bytes()).hexdigest()
                result = Path(sorted_bam(source, None, None, 1))
                self.assertEqual(result, parent / 'input.asmsv.qname.bam')
                receipt = json.loads(Path(str(result) + '.json').read_text())
                self.assertEqual(receipt['source'], identity(source))
                self.assertEqual(receipt['cache'], identity(result))
                self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), before)
                with patch('asmsv_next.cli.subprocess.run', side_effect=AssertionError('cache miss')):
                    self.assertEqual(sorted_bam(source, None, None, 1), str(result))

    def test_stale_missing_and_corrupt_caches_are_rebuilt(self):
        for damage in ('source_changed', 'invalid_json', 'missing_receipt', 'truncated_bam', 'wrong_header'):
            with self.subTest(damage=damage), tempfile.TemporaryDirectory() as root:
                directory = Path(root)
                source = self.make_bam(directory)
                result = Path(sorted_bam(source, None, None, 1))
                receipt = Path(str(result) + '.json')
                if damage == 'source_changed':
                    stat = source.stat()
                    os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
                elif damage == 'invalid_json':
                    receipt.write_text('{interrupted write')
                elif damage == 'missing_receipt':
                    receipt.unlink()
                else:
                    if damage == 'truncated_bam':
                        result.write_bytes(result.read_bytes()[:-28])
                    else:
                        result.write_bytes(source.read_bytes())
                    # Keep metadata consistent to exercise BAM integrity/header validation.
                    saved = json.loads(receipt.read_text()); saved['cache'] = identity(result)
                    receipt.write_text(json.dumps(saved))
                with patch('asmsv_next.cli.subprocess.run', wraps=subprocess.run) as sorting:
                    self.assertEqual(sorted_bam(source, None, None, 1), str(result))
                    self.assertEqual(sorting.call_count, 1)
                pysam.quickcheck(str(result))
                with pysam.AlignmentFile(str(result)) as stream:
                    self.assertEqual([row.query_name for row in stream], ['q1', 'q2', 'q2'])

    def test_explicit_sorted_bam_takes_priority_without_creating_cache(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            source = self.make_bam(directory)
            supplied = sorted_bam(source, None, None, 1)
            work = directory / 'unused'
            with patch('asmsv_next.cli.subprocess.run', side_effect=AssertionError('unnecessary sort')):
                self.assertEqual(sorted_bam(source, supplied, work, 1), supplied)
            self.assertFalse(work.exists())

    def test_explicit_work_directory_reuses_legacy_receipt(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            source = self.make_bam(directory)
            work = directory / 'custom'
            result = Path(sorted_bam(source, None, work, 1))
            self.assertEqual(result.parent, work)
            self.assertFalse((directory / 'input.asmsv.qname.bam').exists())
            receipt = next(work.glob('*.json'))
            receipt.write_text(json.dumps(identity(source)))
            with patch('asmsv_next.cli.subprocess.run', side_effect=AssertionError('legacy cache miss')):
                self.assertEqual(sorted_bam(source, None, work, 1), str(result))

    def test_source_change_during_sort_is_not_published(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            source = self.make_bam(directory)
            original = subprocess.run
            def changed(argv, **kwargs):
                result = original(argv, **kwargs)
                stat = source.stat()
                os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
                return result
            with patch('asmsv_next.cli.subprocess.run', side_effect=changed):
                with self.assertRaisesRegex(ValueError, 'Input BAM changed'):
                    sorted_bam(source, None, None, 1)
            self.assertEqual(list(directory.iterdir()), [source])

    def test_interrupt_removes_sort_chunks_without_publishing_cache(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            work = directory / 'cache'; work.mkdir()
            scratch = directory / 'scratch'; scratch.mkdir()
            source = self.make_bam(directory)
            def interrupted(argv, **kwargs):
                Path(argv[argv.index('-T') + 1] + '.0000.bam').write_bytes(b'partial')
                Path(argv[argv.index('-o') + 1]).write_bytes(b'partial')
                raise KeyboardInterrupt()
            with patch('asmsv_next.cli.subprocess.run', side_effect=interrupted):
                with self.assertRaises(KeyboardInterrupt):
                    sorted_bam(source, None, work, 1, scratch)
            self.assertEqual(list(work.iterdir()), [])
            self.assertEqual(list(scratch.iterdir()), [])

    def test_failed_stage_retains_diagnosis_and_elapsed_time(self):
        with tempfile.TemporaryDirectory() as root:
            record = Path(root) / 'run.json'
            manifest = {'status': 'RUNNING'}
            with self.assertLogs('asmsv_next.cli', level='INFO') as logged:
                with self.assertRaises(KeyboardInterrupt):
                    with run_stage(manifest, record, 'prepare_hap2'):
                        self.assertEqual(json.loads(record.read_text())['stage'], 'prepare_hap2')
                        raise KeyboardInterrupt()
            saved = json.loads(record.read_text())
            self.assertGreaterEqual(saved['stage_seconds']['prepare_hap2'], 0)
            self.assertTrue(any('Stopped prepare_hap2' in message for message in logged.output))


if __name__ == '__main__':
    unittest.main()
