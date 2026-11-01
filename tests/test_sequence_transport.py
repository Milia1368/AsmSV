from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from multiprocessing import shared_memory
import pickle
from types import SimpleNamespace
from unittest.mock import patch
import unittest

from asmsv_next.input_io import QuerySequence
from asmsv_next.sequence_transport import SequenceTransport
from asmsv_next import pipeline
from asmsv_next.config import SVConfig


def child_slices(sequence):
    try:return sequence[37:213],sequence[-33:],sequence[0]
    finally:sequence.close_shared()


def failing_contig(name,records,sequence,length,haplotype,trace_dir):
    try:
        sequence[0]
        raise RuntimeError('worker failure for cleanup test')
    finally:sequence.close_shared()


class SequenceTransportTests(unittest.TestCase):
    def test_shared_slices_spawn_pickle_and_cleanup(self):
        text='aCGTNRYK'*200000
        original=QuerySequence('q',len(text),sequence=text)
        transport=SequenceTransport(minimum=1)
        try:
            shared=transport.prepare(original)
            if shared.shared_name is None:self.skipTest('No shared-memory allocation available')
            self.assertIsNone(shared.sequence)
            self.assertLess(len(pickle.dumps(shared)),2048)
            self.assertEqual(shared[-15:],text[-15:])
            shared.close_shared()
            with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context('spawn')) as pool:
                self.assertEqual(pool.submit(child_slices,shared).result(),(text[37:213],text[-33:],text[0]))
            name=shared.shared_name;transport.release(name)
            with self.assertRaises(FileNotFoundError):shared_memory.SharedMemory(name=name)
        finally:transport.close()

    def test_insufficient_space_and_forced_copy_preserve_original_sequence(self):
        original=QuerySequence('q',5000,sequence='ACGTN'*1000)
        with patch('asmsv_next.sequence_transport.shutil.disk_usage',return_value=SimpleNamespace(free=0)):
            transport=SequenceTransport(minimum=1)
            self.assertIs(transport.prepare(original),original)
            self.assertEqual(transport.copied_bases,5000)
        copied=SequenceTransport('copy',minimum=1)
        self.assertIs(copied.prepare(original),original)
        self.assertFalse(copied.owners)

    def test_allocation_failure_falls_back_to_copy(self):
        original=QuerySequence('q',5000,sequence='ACGTN'*1000)
        with patch('asmsv_next.sequence_transport.shutil.disk_usage',return_value=SimpleNamespace(free=10**12)), \
             patch('asmsv_next.sequence_transport.shared_memory.SharedMemory',side_effect=OSError('no space')):
            transport=SequenceTransport(minimum=1)
            self.assertIs(transport.prepare(original),original)
            self.assertFalse(transport.owners)

    def test_worker_failure_releases_owned_shared_buffers(self):
        original=QuerySequence('q',1600000,sequence='ACGT'*400000)
        transport=SequenceTransport(minimum=1)
        probe=transport.prepare(original)
        if probe.shared_name is None:self.skipTest('No shared-memory allocation available')
        created=[probe.shared_name]
        transport.close()
        def groups(path,label,query_fasta=None,sequence_transform=None):
            sequence=sequence_transform(original);created.append(sequence.shared_name)
            yield 'q',[],sequence,len(sequence),label
        SVConfig.SEQUENCE_TRANSPORT='auto'
        with patch.object(pipeline,'bam_groups',groups),patch.object(pipeline,'call_contig',failing_contig):
            with self.assertRaisesRegex(RuntimeError,'worker failure'):
                pipeline.run_haplotype('unused','hap1',1)
        for name in created:
            if name:
                with self.assertRaises(FileNotFoundError):shared_memory.SharedMemory(name=name)


if __name__=='__main__':unittest.main()
