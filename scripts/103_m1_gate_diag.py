#!/usr/bin/env python3
"""Why does M1 fail outside a species' modal size band? Two stratifications, then stop.

M1's aligned rate is NON-MONOTONIC in insert length -- measured across four
species it peaks at 1200-1599 bp (81-87%), dips at 1600-2399 (20-54%), partly
recovers at 2400-3499, and collapses at >=3500 (0.8-21.7%). A coverage threshold
cannot make that shape. Two candidate causes, and they imply opposite fixes:

  COVERAGE TOO STRICT      identity high, coverage low  -> relax --min-cov
  BOUNDARY DISAGREEMENT    identity high, coverage low  -> switch to CONTAINMENT
                           (is the insert inside an Arm A element?) rather than
                           two-way coverage; relaxing the threshold would only
                           admit mismatches
  SOMETHING ELSE           identity itself low

Arm A derives element boundaries from self-alignment repeat blocks; Arm B
derives the insert from a bubble length difference. Those agree on a species'
dominant families and can disagree systematically elsewhere -- which is a
boundary problem, not a threshold problem.

TEST 2 is free: are the failing short S. aureus inserts the repeat-context ones?
Stratify <800 bp by junction overlap > 15.

    103_m1_gate_diag.py --recon <dirs> --arma <dir> --mobility <tsv> --out x
"""
import argparse, collections, csv, glob, os, statistics, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


def read_alleles(path):
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--arma", required=True)
    ap.add_argument("--mobility", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", type=int, default=400)
    args = ap.parse_args()
    recon = []
    for d in args.recon:
        recon.extend(sorted(glob.glob(d)) if any(c in d for c in "*?[") else [d])

    mob = {r["locus_id"]: r for r in csv.DictReader(open(args.mobility), delimiter="\t")}
    loci, aseq, carr = {}, {}, collections.defaultdict(dict)
    for d in recon:
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            loci[r["locus_id"]] = r
        aseq.update(read_alleles(os.path.join(d, "alleles.fna")))
        ap_ = os.path.join(d, "alleles.tsv")
        if os.path.exists(ap_):
            for r in csv.DictReader(open(ap_), delimiter="\t"):
                carr[r["locus_id"]][r["allele_id"]] = \
                    [g for g in r.get("genomes", "").split(",") if g]

    # ---- TEST 2 (free): short inserts, split by repeat context -------------
    log("=" * 74)
    log("TEST 2 -- are the failing SHORT inserts the repeat-context ones?")
    log("")
    log("  %-22s %8s %10s" % ("<800 bp inserts", "n", "M1 aligned"))
    log("  " + "-" * 44)
    for lab, want in (("overlap <= 15", False), ("overlap > 15  (repeat)", True)):
        v = [r for r in mob.values()
             if int(r["inserted_len"]) < 800
             and (int(loci.get(r["locus_id"], {}).get("junction_ambiguity_bp") or 0) > 15) == want]
        if not v:
            continue
        a = sum(1 for r in v if r["M1_within_genome"] == "True")
        log("  %-22s %8d %9.1f%%" % (lab, len(v), 100.0 * a / len(v)))

    # ---- TEST 1: identity vs coverage on the failing bands -----------------
    log("")
    log("=" * 74)
    log("TEST 1 -- for inserts that FAILED M1, what did the best Arm A hit score?")
    log("")
    bands = [("1200-1599 (modal)", 1200, 1600), ("1600-2399", 1600, 2400),
             (">=3500", 3500, 10 ** 9)]
    tmp = tempfile.mkdtemp(prefix="m1diag.")
    log("  %-20s %6s %9s %9s %9s %9s" %
        ("band", "n", "no hit", "id>=.95", "cov>=.90", "med cov"))
    log("  " + "-" * 68)
    for lab, lo, hi in bands:
        fails = [r for r in mob.values()
                 if lo <= int(r["inserted_len"]) < hi and r["M1_within_genome"] != "True"]
        if not fails:
            continue
        import random
        random.seed(5)
        sel = random.sample(fails, min(args.sample, len(fails)))
        bygen = collections.defaultdict(list)
        for r in sel:
            lid = r["locus_id"]
            row = loci.get(lid)
            if not row:
                continue
            g = (carr.get(lid, {}).get(row["longest_allele"]) or [None])[0]
            if g:
                bygen[g].append(lid)
        best = {}
        for g, lids in bygen.items():
            el = os.path.join(args.arma, g + "_elements.fna")
            if not os.path.exists(el):
                continue
            q = os.path.join(tmp, g + ".fa")
            with open(q, "w") as fh:
                for lid in lids:
                    row = loci[lid]
                    ls = aseq.get(lid, {}).get(row["longest_allele"], "")
                    off, il = int(row["lcp_bp"]), int(row["inserted_len"])
                    if ls:
                        fh.write(">%s\n%s\n" % (lid, ls[off:off + il]))
            paf = q + ".paf"
            with open(paf, "w") as fh:
                subprocess.run(["minimap2", "-c", "-x", "asm10", "-t", "16", el, q],
                               stdout=fh, stderr=subprocess.DEVNULL, check=False)
            for line in open(paf):
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                qn, ql = f[0], int(f[1])
                nm, bl = int(f[9]), int(f[10])
                ident, cov = (nm / bl if bl else 0), (bl / ql if ql else 0)
                if qn not in best or ident * cov > best[qn][0] * best[qn][1]:
                    best[qn] = (ident, cov)
        n = len(sel)
        nohit = sum(1 for r in sel if r["locus_id"] not in best)
        hi_id = sum(1 for k, v in best.items() if v[0] >= 0.95)
        hi_cov = sum(1 for k, v in best.items() if v[1] >= 0.90)
        med = statistics.median([v[1] for v in best.values()]) if best else float("nan")
        log("  %-20s %6d %8.1f%% %8.1f%% %8.1f%% %9.3f"
            % (lab, n, 100.0 * nohit / n, 100.0 * hi_id / n, 100.0 * hi_cov / n, med))
    log("")
    log("  READING: high identity + low coverage = BOUNDARY DISAGREEMENT ->")
    log("  switch M1 to CONTAINMENT (insert inside an element), not two-way")
    log("  coverage. Low identity = the insert is simply not that element, and")
    log("  relaxing the gate would admit mismatches.")
    log("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
