"""Bounded shared sequence transport; callers see the same string slices.

The parent owns/unlinks buffers after their futures complete. Workers only map
them and close their handles. Insufficient shared-memory space uses the original
copy transport, without changing any sequence or variant rule.
"""
from dataclasses import replace
from multiprocessing import shared_memory
from pathlib import Path
import os
import shutil


class SequenceTransport:
    def __init__(self, mode='auto', minimum=1024*1024):
        self.mode=mode;self.minimum=minimum;self.owners={}
        self.shared_bases=0;self.copied_bases=0;self.shared_queries=0

    def prepare(self, sequence):
        if sequence.sequence is None:return sequence
        n=len(sequence)
        if self.mode=='copy' or n<self.minimum:
            self.copied_bases+=n;return sequence
        owner=None
        try:
            # /dev/shm is bounded on Linux; leave room for other processes and
            # avoid a successful ftruncate followed by SIGBUS while copying.
            if not Path('/dev/shm').is_dir() or shutil.disk_usage('/dev/shm').free < n+n//4+32*1024*1024:
                self.copied_bases+=n;return sequence
            owner=shared_memory.SharedMemory(create=True,size=n)
            if hasattr(os,'posix_fallocate'):
                fd=os.open(str(Path('/dev/shm')/owner.name),os.O_RDWR)
                try:os.posix_fallocate(fd,0,n)
                finally:os.close(fd)
        except OSError:
            if owner is not None:owner.close();owner.unlink()
            self.copied_bases+=n;return sequence
        try:
            # Bounded temporary encoding, rather than another chromosome-sized
            # bytes object in addition to the original string.
            step=1024*1024
            for start in range(0,n,step):
                chunk=sequence.sequence[start:start+step].encode('ascii')
                owner.buf[start:start+len(chunk)]=chunk
            self.owners[owner.name]=owner
            self.shared_bases+=n;self.shared_queries+=1
            return replace(sequence,sequence=None,shared_name=owner.name)
        except BaseException:
            owner.close();owner.unlink();raise

    def release(self, name):
        owner=self.owners.pop(name,None)
        if owner is not None:
            owner.close();owner.unlink()

    def close(self):
        for name in list(self.owners):self.release(name)
