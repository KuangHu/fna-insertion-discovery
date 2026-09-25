"""Minimal PAF reader. Stdlib only."""
from dataclasses import dataclass


@dataclass(slots=True)
class Paf:
    qname: str
    qlen: int
    qs: int
    qe: int
    strand: str
    tname: str
    tlen: int
    ts: int
    te: int
    nmatch: int
    alen: int
    mapq: int
    tags: dict

    @property
    def identity(self) -> float:
        """Gap-compressed identity if de/dv tag present, else blast identity."""
        if "de" in self.tags:
            return (1.0 - float(self.tags["de"])) * 100.0
        if "dv" in self.tags:
            return (1.0 - float(self.tags["dv"])) * 100.0
        return 100.0 * self.nmatch / self.alen if self.alen else 0.0

    @property
    def qspan(self) -> int:
        return self.qe - self.qs

    @property
    def tspan(self) -> int:
        return self.te - self.ts


def _tag(tok: str):
    name, typ, val = tok.split(":", 2)
    if typ == "i":
        return name, int(val)
    if typ == "f":
        return name, float(val)
    return name, val


def parse(path):
    """Yield Paf records."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tags = {}
            for tok in f[12:]:
                if tok.count(":") >= 2:
                    k, v = _tag(tok)
                    tags[k] = v
            yield Paf(f[0], int(f[1]), int(f[2]), int(f[3]), f[4],
                      f[5], int(f[6]), int(f[7]), int(f[8]),
                      int(f[9]), int(f[10]), int(f[11]), tags)
