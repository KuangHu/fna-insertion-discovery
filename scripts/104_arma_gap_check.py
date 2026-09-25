#!/usr/bin/env python3
"""Is M1's no-hit a TRUE negative, or an Arm A family-construction gap?

76-91% of M1 failures produce no Arm A hit at all -- not a low-coverage hit, a
hit of any kind. That has two readings with opposite consequences:

  TRUE NEGATIVE     the insert really is single-copy in that genome; M1 is right
  ARM A GAP         it IS multi-copy, but Arm A never grouped it into a family

The modal size band forces the question: 90.8% no-hit sits exactly where each
species' dominant, high-copy elements live, and "dominant element that is not
multi-copy in its own genome" is close to a contradiction.

So bypass Arm A entirely: map each no-hit insert straight back against its own
carrier genome and COUNT THE COPIES. Arm A's own criterion is >=3 copies at
distinct loci (--min-copies 3, --min-locus-sep), so the same bar is applied here
with no family model in between.

    < 3 copies  -> true negative, M1 is correct, nothing to fix
   >= 3 copies  -> Arm A missed a family; the M line's ONLY source has a gap,
                   which is upstream of every catalogue produced so far

    104_arma_gap_check.py --recon <dirs> --mobility <tsv> --genome-dir <dir> \\
        --lo 1200 --hi 1600 --out x
"""
import argparse, collections, csv, glob, os, random, subprocess, sys, tempfile
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
    ap.add_argument("--mobility", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--lo", type=int, default=1200)
    ap.add_argument("--hi", type=int, default=1600)
    ap.add_argument("--sample", type=int, default=300)
    ap.add_argument("--min-ident", type=float, default=0.95)
    ap.add_argument("--min-cov", type=float, default=0.90)
    ap.add_argument("--min-sep", type=int, default=1000)
    ap.add_argument("--out", required=True)
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
        a = os.path.join(d, "alleles.tsv")
        if os.path.exists(a):
            for r in csv.DictReader(open(a), delimiter="\t"):
                carr[r["locus_id"]][r["allele_id"]] = \
                    [g for g in r.get("genomes", "").split(",") if g]

    fails = [r for r in mob.values()
             if args.lo <= int(r["inserted_len"]) < args.hi
             and r["M1_within_genome"] != "True"]
    random.seed(9)
    sel = random.sample(fails, min(args.sample, len(fails)))
    log("%d M1 failures in %d-%d bp; sampling %d" % (len(fails), args.lo, args.hi, len(sel)))

    bygen = collections.defaultdict(list)
    for r in sel:
        lid = r["locus_id"]
        row = loci.get(lid)
        if not row:
            continue
        g = (carr.get(lid, {}).get(row["longest_allele"]) or [None])[0]
        if g:
            bygen[g].append(lid)
    log("across %d carrier genomes" % len(bygen))

    tmp = tempfile.mkdtemp(prefix="armagap.")
    counts = {}
    for g, lids in bygen.items():
        fa = os.path.join(args.genome_dir, g + ".fna")
        if not os.path.exists(fa):
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
            subprocess.run(["minimap2", "-c", "-x", "asm10", "-N", "60", "-p", "0.1",
                            "-t", "16", fa, q], stdout=fh,
                           stderr=subprocess.DEVNULL, check=False)
        hit = collections.defaultdict(list)
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qn, ql = f[0], int(f[1])
            nm, bl = int(f[9]), int(f[10])
            if bl == 0 or ql == 0:
                continue
            if nm / bl >= args.min_ident and bl / ql >= args.min_cov:
                hit[qn].append((f[5], int(f[7])))
        for lid in lids:
            # distinct loci: same contig within --min-sep counts once
            seen = []
            for c, s in sorted(hit.get(lid, [])):
                if not any(c == pc and abs(s - ps) < args.min_sep for pc, ps in seen):
                    seen.append((c, s))
            counts[lid] = len(seen)

    h = collections.Counter()
    for lid in sel:
        n = counts.get(lid["locus_id"], 0) if isinstance(lid, dict) else 0
        h["0-1" if n <= 1 else "2" if n == 2 else ">=3"] += 1
    tot = sum(h.values())
    log("=" * 74)
    log("ARM A GAP CHECK -- copies counted DIRECTLY in the carrier genome")
    log("  (Arm A's own bar: >=3 copies at distinct loci, >=%d bp apart)" % args.min_sep)
    log("")
    for k in ("0-1", "2", ">=3"):
        log("  %-6s copies  %6d  %5.1f%%" % (k, h[k], 100.0 * h[k] / tot if tot else 0))
    log("")
    if tot and h[">=3"] / tot > 0.25:
        log("  *** ARM A GAP: %.1f%% of these ARE multi-copy but were never" % (100.0 * h[">=3"] / tot))
        log("  grouped into a family. The M line's only source under-reports,")
        log("  and that is upstream of every catalogue produced so far. ***")
    else:
        log("  TRUE NEGATIVES: these inserts really are single- or low-copy in")
        log("  their own genome. M1 is correct and the gate needs no change.")
    log("=" * 74)
    with open(args.out + "_armagap.tsv", "w") as fh:
        fh.write("locus_id\tcopies_direct\n")
        for lid, n in sorted(counts.items()):
            fh.write("%s\t%d\n" % (lid, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
