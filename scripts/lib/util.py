"""Union-find, interval merging, small shared helpers. Stdlib only."""
import hashlib
import subprocess
import sys
import time


class UnionFind:
    def __init__(self):
        self.p = {}

    def add(self, x):
        self.p.setdefault(x, x)

    def find(self, x):
        self.add(x)
        r = x
        while self.p[r] != r:
            r = self.p[r]
        while self.p[x] != r:
            self.p[x], x = r, self.p[x]
        return r

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra

    def groups(self):
        g = {}
        for x in list(self.p):
            g.setdefault(self.find(x), []).append(x)
        return list(g.values())


def merge_intervals(ivs, slack=0):
    """ivs: list of (start, end, payload...). Returns merged (start, end)."""
    if not ivs:
        return []
    s = sorted((i[0], i[1]) for i in ivs)
    out = [list(s[0])]
    for a, b in s[1:]:
        if a <= out[-1][1] + slack:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def jaccard(a0, a1, b0, b1):
    ov = overlap(a0, a1, b0, b1)
    if ov == 0:
        return 0.0
    union = (a1 - a0) + (b1 - b0) - ov
    return ov / union if union else 0.0


def seq_md5(seq):
    return hashlib.md5(seq.upper().encode()).hexdigest()


def run(cmd, **kw):
    """Run a command, raising with stderr on failure."""
    kw.setdefault("check", True)
    kw.setdefault("capture_output", True)
    kw.setdefault("text", True)
    try:
        return subprocess.run(cmd, **kw)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("FAILED: %s\n%s\n" % (" ".join(map(str, cmd)), e.stderr))
        raise


def log(msg):
    sys.stderr.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))
    sys.stderr.flush()
