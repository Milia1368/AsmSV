"""Content-addressed insertion FASTA, published only after successful VCF writing."""
import hashlib
import os
from pathlib import Path
import uuid
import pysam


class SequenceStore:
    def __init__(self,path):
        self.path=Path(path);self.temporary=self.path.with_name(self.path.name+'.'+uuid.uuid4().hex+'.tmp')
        self.stream=self.temporary.open('w');self.seen={};self.finished=False

    def annotate(self,candidate):
        if candidate.type not in ('INS','DELINS'):return
        sequence=candidate.info.get('SEQ','')
        if not sequence:return
        digest=hashlib.sha256(sequence.encode()).hexdigest();identifier='ASMSVSEQ.'+digest
        if digest not in self.seen:
            self.stream.write('>'+identifier+' length='+str(len(sequence))+' sha256='+digest+'\n')
            for start in range(0,len(sequence),80):self.stream.write(sequence[start:start+80]+'\n')
            self.seen[digest]=len(sequence)
        elif self.seen[digest]!=len(sequence):
            raise ValueError('Insertion sequence digest collision')
        candidate.info.update(SEQID=identifier,SEQLEN=len(sequence),SEQSHA256=digest,SEQFILE=self.path.name)

    def finish(self):
        if self.finished:return
        self.stream.close();os.replace(self.temporary,self.path)
        index=Path(str(self.path)+'.fai')
        if self.seen:pysam.faidx(str(self.path))
        elif index.exists():index.unlink()
        self.finished=True

    def abort(self):
        if not self.stream.closed:self.stream.close()
        if self.temporary.exists():self.temporary.unlink()
