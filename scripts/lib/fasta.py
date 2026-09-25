"""Random-access FASTA via a samtools-style .fai index. Stdlib only."""
import os
import subprocess

_RC = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.encode().translate(_RC)[::-1].decode()


class Fasta:
    """Coordinate convention: 0-based, half-open [start, end)."""

    def __init__(self, path, samtools="samtools"):
        self.path = str(path)
        fai = self.path + ".fai"
        if not os.path.exists(fai):
            subprocess.run([samtools, "faidx", self.path], check=True)
        self.idx = {}
        with open(fai) as fh:
            for line in fh:
                name, ln, off, lb, wb = line.split("\t")[:5]
                self.idx[name] = (int(ln), int(off), int(lb), int(wb))
        self._fh = open(self.path, "rb")

    def __contains__(self, name):
        return name in self.idx

    def names(self):
        return list(self.idx)

    def length(self, name):
        return self.idx[name][0]

    def total(self):
        return sum(v[0] for v in self.idx.values())

    def fetch(self, name, start=0, end=None):
        ln, off, lb, wb = self.idx[name]
        if end is None or end > ln:
            end = ln
        start = max(0, start)
        if start >= end:
            return ""
        # byte offset of the base at `start`
        self._fh.seek(off + start // lb * wb + start % lb)
        need = end - start
        # read enough bytes to cover newlines
        raw = self._fh.read(need + need // lb + 2)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:need]
        return seq.decode()

    def close(self):
        self._fh.close()


def write_fasta(fh, name, seq, width=80):
    fh.write(">" + name + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")
