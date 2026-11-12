"""Random-access FASTA via a samtools-style .fai index. Stdlib only."""
import os
import subprocess

# full IUPAC alphabet: R/Y, K/M, B/V, D/H swap; S, W, N are self-complementary
_RC = bytes.maketrans(b"ACGTRYKMBVDHSWNacgtrykmbvdhswn",
                      b"TGCAYRMKVBHDSWNtgcayrmkvbhdswn")
_NOT_BASES = (b">", b" ", b"\t")     # never inside a sequence line


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
            for n, line in enumerate(fh, 1):
                f = line.rstrip("\n").split("\t")
                try:
                    name, ln, off, lb, wb = f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4])
                except (IndexError, ValueError):
                    raise ValueError("%s line %d: malformed .fai entry" % (fai, n))
                if ln < 0 or off < 0 or lb < 1 or wb < lb or name in self.idx:
                    raise ValueError("%s line %d: invalid .fai entry %r" % (fai, n, line[:80]))
                self.idx[name] = (ln, off, lb, wb)
        self._fh = open(self.path, "rb")
        self._verify_index(fai)

    def _verify_index(self, fai):
        """Bind the .fai to THIS file, record by record in file order.

        Walking the entries by offset, each record's header line must start
        exactly where the previous record ended (blank lines allowed), name the
        entry, and end at the entry's offset. Then EVERY line of the record is
        checked against the index (_check_layout): each line but the last holds
        exactly line_bases bases followed by one line end of line_bytes -
        line_bases bytes, and the last line holds the remaining bases. After the
        last record only line ends may remain. So an empty, truncated or stale
        index -- an unindexed record, a sequence longer or shorter than indexed,
        a header that moved, a line of the wrong length anywhere in a record --
        fails here instead of fetch() silently returning bases from the wrong
        place.

        Every open reads the whole file once (~12 ms for a 5 MB genome in the
        page cache, measured 2026-11-12). Nothing is cached across opens: file
        timestamps on /global/scratch have 1 s resolution, so no cheap key can
        tell an in-place rewrite from the file that was checked."""
        size = os.path.getsize(self.path)
        if not self.idx:
            if size:
                raise ValueError("empty .fai for non-empty %s" % self.path)
            return
        prev_end = 0                 # where the next record (or EOF) may begin
        for name, (ln, off, lb, wb) in sorted(self.idx.items(), key=lambda kv: kv[1][1]):
            bad = "stale or wrong .fai for %s: entry %r" % (self.path, name)
            if off < prev_end:
                raise ValueError(bad + " (overlaps the previous record)")
            self._fh.seek(prev_end)
            gap = self._fh.read(off - prev_end).lstrip(b"\r\n")
            if (not gap.startswith(b">") or not gap.endswith(b"\n")
                    or gap.count(b"\n") != 1):
                raise ValueError(bad + " (offset does not follow exactly one header line; "
                                 "an unindexed record or a moved header?)")
            if gap[1:].rstrip(b"\r\n").split()[:1] != [name.encode()]:
                raise ValueError(bad + " (header line names a different record)")
            prev_end = self._check_layout(bad, ln, off, lb, wb, size)
        tail_len = size - prev_end
        self._fh.seek(prev_end)
        if tail_len > 4096 or self._fh.read(tail_len).strip(b"\r\n"):
            raise ValueError("stale or wrong .fai for %s: the file continues after the last "
                             "indexed base (unindexed record or longer sequence)" % self.path)

    def _check_layout(self, bad, ln, off, lb, wb, size):
        """Every line of one record, read in chunks of whole lines (bounded
        memory). Returns the offset just past the record's last base.

        A full line is lb bases then a line end: '\\n' when wb - lb == 1,
        '\\r\\n' when wb - lb == 2 (the only two samtools supports). The line
        ends must sit at exactly the positions lb, lb + wb, lb + 2*wb, ... and
        nowhere else, and no '>' or blank may appear among the bases. The last
        line holds the remaining 1..lb bases and no line end; what follows it
        (the next header, or EOF) is checked by the caller."""
        if ln == 0:
            return off
        nfull = (ln - 1) // lb       # lines that must hold exactly lb bases
        rest = ln - nfull * lb       # bases on the last line, 1..lb
        end = off + nfull * wb + rest
        if end > size:
            raise ValueError(bad + " (sequence runs past end of file)")
        eol = {1: b"\n", 2: b"\r\n"}.get(wb - lb)
        if nfull and eol is None:
            raise ValueError(bad + " (line_bytes - line_bases must be 1 or 2, got %d)"
                             % (wb - lb))
        self._fh.seek(off)
        per = max(1, (1 << 22) // wb)            # ~4 MB of whole lines per read
        done = 0
        while done < nfull:
            k = min(per, nfull - done)
            block = self._fh.read(k * wb)
            ok = (len(block) == k * wb
                  and block.count(b"\n") == k
                  and block.count(b"\r") == (k if len(eol) == 2 else 0)
                  and all(block[lb + i::wb] == eol[i:i + 1] * k for i in range(len(eol)))
                  and not any(c in block for c in _NOT_BASES))
            if not ok:
                raise ValueError(bad + " (sequence lines %d-%d do not match line_bases=%d, "
                                 "line_bytes=%d)" % (done + 1, done + k, lb, wb))
            done += k
        last = self._fh.read(rest)
        if (len(last) != rest or b"\n" in last or b"\r" in last
                or any(c in last for c in _NOT_BASES)):
            raise ValueError(bad + " (last line does not hold the %d remaining bases)" % rest)
        return end

    def __contains__(self, name):
        return name in self.idx

    def names(self):
        return list(self.idx)

    def length(self, name):
        return self.idx[name][0]

    def total(self):
        return sum(v[0] for v in self.idx.values())

    def fetch(self, name, start=0, end=None):
        """Bases [start, end) of `name`, clamped to the contig.

        The byte range is computed exactly from the index (line_bases,
        line_bytes), so CRLF files (line_bytes = line_bases + 2) read correctly,
        and the result must have exactly end - start bases: a truncated file or
        a stale .fai raises instead of returning a shorter sequence."""
        ln, off, lb, wb = self.idx[name]
        if end is None or end > ln:
            end = ln
        start = max(0, start)
        if start >= end:
            return ""
        b0 = off + start // lb * wb + start % lb
        b1 = off + (end - 1) // lb * wb + (end - 1) % lb + 1
        self._fh.seek(b0)
        raw = self._fh.read(b1 - b0)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"")
        if len(seq) != end - start:
            raise IOError("short read from %s: %s:%d-%d gave %d of %d bases "
                          "(truncated file or stale .fai?)"
                          % (self.path, name, start, end, len(seq), end - start))
        return seq.decode()

    def close(self):
        self._fh.close()


def write_fasta(fh, name, seq, width=80):
    fh.write(">" + name + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")
