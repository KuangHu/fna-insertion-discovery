#!/usr/bin/env python3
"""AUDIT A, redone at OCCURRENCE grain against the reference each row actually used.

The first pass asked a weaker question in a weaker place.

  WEAKER PLACE. It walked the locations table with `GROUP BY db_id`, which
  keeps ONE arbitrary coord_frame per event. An event may be supported by
  several empty references, and the one SQLite happened to return need not be
  the one any official occurrence used. An anomaly in reference X does not
  impeach an occurrence whose empty site was taken from reference Y. This pass
  keys on (empty_accession, empty_contig, empty_junction_L) AS RECORDED ON THE
  OCCURRENCE ROW, so every alarm lands on a row that depends on it.

  WEAKER QUESTION. It matched a ONE-SIDED probe, host-left + insert-start, and
  reported the distance to the recorded junction. A one-sided 40+40 match shows
  only that the insert's first bases occur after that host context somewhere;
  repeats and shared element termini produce exactly that. It cannot
  distinguish "this assembly carries the element here" from "this 40-mer is
  ambiguous". This pass tests the STRUCTURE: left junction, the full span the
  element would occupy, and the right end, in one orientation at a time.

A hit far from the recorded junction is NOT an anomaly -- a mobile element is
expected elsewhere in the same genome. Only the recorded locus is examined.

Output is per occurrence so exclusions can be unioned by occurrence_id.
"""
import argparse
import collections
import os
import sys
import zlib

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"
V2 = W + "/inferred_empty_v2/occurrences.tsv"
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

K = 40            # junction k-mer
LEFT_MAX_MM = 2   # mismatches tolerated in the left junction k-mer
ID_MIN = 0.90     # span identity to call the element present
COV_MIN = 0.90    # span coverage to call the element present
REPEAT_N = 5      # genome-wide copies of the left k-mer that make it ambiguous


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    out, k, buf = {}, None, []
    with open(p) as fh:
        for line in fh:
            if line[0] == ">":
                if k is not None:
                    out[k] = "".join(buf).upper()
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if k is not None:
        out[k] = "".join(buf).upper()
    return out


def load_map(p):
    m = {}
    with open(p) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                m[f[0]] = f[1]
    return m


def hamming_le(a, b, cap):
    """True if hamming(a,b) <= cap, with early exit. Lengths must match."""
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > cap:
                return False, n
    return True, n


def identity(a, b):
    """Ungapped identity over the overlapping prefix."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(x == y for x, y in zip(a[:n], b[:n]))
    return same / float(n)


def agree_prefix(sub, ele, floor=0.90):
    """Longest prefix n of `ele` matching `sub` at >= `floor` identity.

    A 40-mer match at the junction followed by random sequence is not an
    element: this says HOW FAR the agreement actually runs before it stops,
    which is the number that separates a shared junction motif from a real
    insertion.
    """
    n = min(len(sub), len(ele))
    mm = 0
    best = 0
    for i in range(n):
        if sub[i] != ele[i]:
            mm += 1
        if (i + 1 - mm) >= floor * (i + 1):
            best = i + 1
    return best


def count_kmer(ctgs, kmer, cap=50):
    """Genome-wide copies of kmer on both strands, capped."""
    n = 0
    rk = revcomp(kmer)
    pats = (kmer,) if rk == kmer else (kmer, rk)
    for v in ctgs.values():
        for p in pats:
            i = v.find(p)
            while i >= 0:
                n += 1
                if n >= cap:
                    return n
                i = v.find(p, i + 1)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--all-rows", action="store_true",
                    help="examine withheld rows too (default: official only)")
    a = ap.parse_args()

    man = load_map(SRC + "/genome_map.tsv")
    ins = read_fasta(SRC + "/insert.fna")

    # ---- group the rows this shard owns by the EMPTY reference they used ----
    by = collections.defaultdict(list)
    n_rows = 0
    with open(V2) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not a.all_rows and f[ix["official_fasta"]] != "yes":
                continue
            acc = f[ix["empty_accession"]]
            if not acc or acc == ".":
                continue
            # Python's str hash is salted per process (PYTHONHASHSEED),
            # so it cannot be used to partition work across processes.
            if zlib.crc32(acc.encode()) % a.nshard != a.shard % a.nshard:
                continue
            jl = f[ix["empty_junction_L"]]
            if jl in (".", ""):
                continue
            by[acc].append((f[ix["occurrence_id"]], f[ix["db_id"]],
                            f[ix["species"]], f[ix["empty_contig"]], int(jl),
                            int(f[ix["inserted_len"]]),
                            f[ix["official_fasta"]]))
            n_rows += 1

    out = open("%s/auditA_occ_%03d.tsv" % (a.out, a.shard), "w")
    out.write("occurrence_id\tdb_id\tspecies\tempty_accession\tempty_contig\t"
              "empty_junction_L\tinserted_len\tofficial\torientation\t"
              "left_mismatches\tagree_prefix_bp\tagree_frac\tspan_identity\t"
              "span_coverage\tright_end_exact\tjunction_kmer_copies\tverdict\n")
    cls = collections.Counter()
    missing = collections.Counter()
    accs = sorted(by)
    for gi, acc in enumerate(accs, 1):
        p = man.get(acc)
        if not p:
            missing["genome_not_in_manifest"] += len(by[acc])
            continue
        try:
            ctgs = read_fasta(p)
        except (IOError, OSError):
            missing["genome_unreadable"] += len(by[acc])
            continue
        for occ, db_id, sp, ctg, jl, ilen, off in by[acc]:
            S = ctgs.get(ctg)
            I = ins.get(db_id, "")
            if S is None:
                missing["contig_not_found"] += 1
                continue
            if not I:
                missing["insert_not_found"] += 1
                continue
            L = len(I)
            if jl < 0 or jl + K > len(S):
                cls["locus_at_contig_end"] += 1
                continue
            win = S[jl:jl + K]
            best = None
            # Orientation is a property of the whole locus: test each frame of
            # the element as a unit, never mix one frame's left with another's.
            for name, E in (("same", I), ("flip", revcomp(I))):
                ok, mm = hamming_le(win, E[:K], LEFT_MAX_MM)
                if ok and (best is None or mm < best[2]):
                    best = (name, E, mm)
            if best is None:
                cls["locus_empty"] += 1
                continue
            name, E, mm = best
            span = S[jl:jl + L]
            cov = len(span) / float(L)
            idt = identity(span, E)
            rexact = (len(span) == L and span[-K:] == E[-K:])
            # Count copies of the GENOME'S OWN junction k-mer, not the
            # element's: with mismatches tolerated the element k-mer may not
            # occur exactly even here, and the question is whether THIS
            # context is repetitive in THIS assembly.
            copies = count_kmer(ctgs, win)
            agree = agree_prefix(span, E)
            if cov < COV_MIN:
                v = "locus_truncated_at_contig_end"
            elif idt >= ID_MIN and rexact:
                v = "element_present_at_locus"
            elif idt >= ID_MIN:
                v = "element_present_right_end_inexact"
            elif agree >= 0.5 * L:
                v = "partial_element_at_locus"
            elif copies >= REPEAT_N:
                v = "left_junction_only_repeat_ambiguous"
            else:
                v = "left_junction_only_divergent_span"
            cls[v] += 1
            out.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d"
                      "\t%.4f\t%.4f\t%.4f\t%s\t%d\t%s\n"
                      % (occ, db_id, sp, acc, ctg, jl, ilen, off, name, mm,
                         agree, agree / float(L), idt, cov,
                         "yes" if rexact else "no", copies, v))
        if gi % 200 == 0:
            sys.stderr.write("  %d/%d genomes\n" % (gi, len(accs)))
    out.close()
    print("shard %d: rows %d  genomes %d" % (a.shard, n_rows, len(accs)))
    for k, v in sorted(cls.items()):
        print("   %-40s %7d" % (k, v))
    for k, v in sorted(missing.items()):
        print("   MISSING %-32s %7d" % (k, v))


if __name__ == "__main__":
    main()
