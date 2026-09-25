#!/usr/bin/env python3
"""Two nulls the shared-target survivors must clear before any claim is made.

TEST 1 -- COINCIDENCE. "Two insertions at the identical base" is only evidence of
target-site preference if coincidence would not produce it. IS elements have real
positional preferences (AT-rich, bent DNA, motifs), so hot sites exist and a few
identical-offset pairs are exactly what preference predicts WITHOUT independent
arrival at one site being remarkable. Two nulls, because they differ in strength:

  NULL A -- uniform: each insertion offset drawn uniformly across the anchored
    interval. Tests "no positional structure at all". Weak, and easy to beat.
  NULL B -- shuffled empirical: offsets resampled from the POOLED observed
    relative-offset distribution, so whatever positional concentration the data
    actually has is preserved, and only the PAIRING within a locus is broken.
    This is the null that matters. If hot sites drive the result, B reproduces it.

TEST 2 -- DIVERGENCE SHAPE. Two elements of near-identical multi-kb length,
83-87% identical, at the identical base is ALSO the signature of one element pair
that diverged in place, or of a composite sharing internal architecture that the
containment test misses because neither middle contains the other. A single
global identity cannot separate these. Windowed identity along the local
alignment can: UNIFORM divergence across kb reads as one ancestral element that
diverged; BLOCKY divergence (alternating conserved and unrelated stretches) reads
as composite.

    97_shared_target_nulls.py --scan shared_target/x_shared_target.tsv \\
        --recon <recon dirs> --out shared_target/x_nulls
"""
import argparse
import collections
import csv
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_alleles_fna(path):
    out, cur = collections.defaultdict(dict), None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur:
                out[cur[0]][cur[1]] = []
        elif cur:
            out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "").upper() for a, v in d.items()}
            for k, d in out.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", required=True, help="_shared_target.tsv from 96")
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--offset-tol", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    random.seed(args.seed)
    rows = list(csv.DictReader(open(args.scan), delimiter="\t"))
    log("%d candidate loci in scan" % len(rows))

    def fnum(r, k):
        try:
            return float(r[k])
        except (ValueError, KeyError):
            return float("nan")

    # interval length = the site's shortest allele; recover it from the alleles
    aseq = {}
    for d in args.recon:
        aseq.update(read_alleles_fna(os.path.join(d, "alleles.fna")))
    L = {}
    for r in rows:
        av = aseq.get(r["locus"], {})
        if av:
            L[r["locus"]] = min(len(s) for s in av.values())

    usable = [r for r in rows if r["locus"] in L and L[r["locus"]] > 2 * args.offset_tol]
    log("%d loci with a usable interval length" % len(usable))

    # ---- TEST 1 ---------------------------------------------------------
    def delta_le(rs, tol):
        return sum(1 for r in rs if abs(int(r["offset1"]) - int(r["offset2"])) <= tol)

    def delta_eq0(rs):
        return sum(1 for r in rs if int(r["offset1"]) == int(r["offset2"]))

    pooled_rel = []
    for r in usable:
        n = L[r["locus"]]
        for k in ("offset1", "offset2"):
            pooled_rel.append(min(1.0, max(0.0, int(r[k]) / n)))

    def permute(rs, mode):
        hit_tol = hit_0 = 0
        for r in rs:
            n = L[r["locus"]]
            if mode == "uniform":
                o1, o2 = random.randrange(n), random.randrange(n)
            else:
                o1 = int(random.choice(pooled_rel) * n)
                o2 = int(random.choice(pooled_rel) * n)
            if abs(o1 - o2) <= args.offset_tol:
                hit_tol += 1
            if o1 == o2:
                hit_0 += 1
        return hit_tol, hit_0

    def run_test(rs, title):
        if not rs:
            return
        obs_t, obs_0 = delta_le(rs, args.offset_tol), delta_eq0(rs)
        log("")
        log("  %s  (n=%d loci)" % (title, len(rs)))
        log("    OBSERVED   delta<=%d : %d      delta==0 : %d"
            % (args.offset_tol, obs_t, obs_0))
        for mode, label in (("uniform", "NULL A uniform"),
                            ("shuffled", "NULL B shuffled-empirical")):
            dt, d0 = [], []
            for _ in range(args.n_perm):
                a, b = permute(rs, mode)
                dt.append(a)
                d0.append(b)
            pt = (sum(1 for x in dt if x >= obs_t) + 1) / (args.n_perm + 1)
            p0 = (sum(1 for x in d0 if x >= obs_0) + 1) / (args.n_perm + 1)
            log("    %-26s delta<=%d mean %5.2f  p=%.4f   |  delta==0 mean %5.2f  p=%.4f"
                % (label, args.offset_tol, statistics.mean(dt), pt,
                   statistics.mean(d0), p0))

    log("=" * 96)
    log("TEST 1 -- COINCIDENCE (%d permutations)" % args.n_perm)
    run_test(usable, "ALL candidates")
    dh = [r for r in usable if r["klass"] == "diverged_homologue"]
    run_test(dh, "diverged homologues")
    surv = [r for r in usable if r["klass"] == "diverged_homologue"
            and r["same_site"] == "True" and r["nested"] == "False"
            and r["well_covered"] == "True"]
    run_test(surv, "SURVIVORS (diverged + well-covered + same-site + not nested)")

    # ---- TEST 2 ---------------------------------------------------------
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                         open_gap_score=-12, extend_gap_score=-0.4)
    log("")
    log("=" * 96)
    log("TEST 2 -- DIVERGENCE SHAPE along the alignment (%d bp windows)" % args.window)
    log("")
    log("  %-20s %7s %7s %8s %8s %8s %8s  %s"
        % ("locus", "len", "global", "win_med", "win_min", "win_max", "CV", "reading"))
    log("  " + "-" * 100)
    out = []
    for r in surv:
        av = aseq.get(r["locus"], {})
        if not av:
            continue
        two = sorted(av.items(), key=lambda kv: -len(kv[1]))[:2]
        s1, s2 = two[0][1], two[1][1]
        n = min(len(s1), len(s2))
        p = 0
        while p < n and s1[p] == s2[p]:
            p += 1
        q = 0
        while q < n - p and s1[len(s1) - 1 - q] == s2[len(s2) - 1 - q]:
            q += 1
        m1, m2 = s1[p:len(s1) - q], s2[p:len(s2) - q]
        if r["orientation"] == "revcomp":
            m2 = rc(m2)
        try:
            a = al.align(m1, m2)[0]
        except Exception:
            continue
        A, B = a.aligned
        # per-column match vector along the aligned blocks only
        col = []
        for (x0, x1), (y0, y1) in zip(A, B):
            col.extend(1 if u == v else 0 for u, v in zip(m1[x0:x1], m2[y0:y1]))
        if len(col) < 3 * args.window:
            continue
        wins = [statistics.mean(col[i:i + args.window])
                for i in range(0, len(col) - args.window + 1, args.window)]
        g = statistics.mean(col)
        cv = statistics.pstdev(wins) / g if g else float("nan")
        lowfrac = sum(1 for w in wins if w < 0.70) / len(wins)
        hifrac = sum(1 for w in wins if w > 0.95) / len(wins)
        # The MEDIAN window, not the blockiness, carries the interpretation.
        # "Blocky" alone conflates two OPPOSITE readings: mostly-identical with
        # one divergent block (same element, internally recombined -- ONE
        # arrival) versus mostly-divergent with one conserved block (two
        # different elements sharing a domain -- consistent with two arrivals).
        wmed = statistics.median(wins)
        if cv < 0.15 and wmed >= 0.90:
            reading = "one element, diverged uniformly - ONE arrival"
        elif wmed >= 0.90:
            reading = "same element + divergent block - ONE arrival, recombined"
        elif wmed <= 0.70:
            reading = "two different elements, shared segment - two arrivals OK"
        else:
            reading = "intermediate - unresolved"
        log("  %-20s %7d %7.3f %8.3f %8.3f %8.3f %8.3f  %s"
            % (r["locus"], len(col), g, statistics.median(wins), min(wins),
               max(wins), cv, reading))
        out.append((r["locus"], len(col), g, cv, lowfrac, hifrac, reading, wins))
    log("")
    log("  READING KEY -- the MEDIAN window identity decides, not the blockiness:")
    log("    median >= 0.90  the two middles are THE SAME ELEMENT over most of")
    log("                    their length. Identical offset then needs no")
    log("                    coincidence at all: one ancestral insertion,")
    log("                    later internal recombination. NOT two arrivals,")
    log("                    and the global identity was hiding it.")
    log("    median <= 0.70  mostly divergent with a shared segment -- two")
    log("                    different elements sharing a domain. Consistent")
    log("                    with two arrivals at one site.")
    with open(args.out + "_divergence_shape.tsv", "w") as fh:
        fh.write("locus\taln_cols\tglobal_identity\twindow_CV\tfrac_win_lt0.70"
                 "\tfrac_win_gt0.95\treading\twindow_identities\n")
        for lo, nc, g, cv, lf, hf, rd, wins in out:
            fh.write("%s\t%d\t%.4f\t%.4f\t%.3f\t%.3f\t%s\t%s\n"
                     % (lo, nc, g, cv, lf, hf, rd,
                        ",".join("%.3f" % w for w in wins)))
    log("=" * 96)
    log("wrote %s_divergence_shape.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
