#!/usr/bin/env python3
"""Within-genome copy structure: four fields, not one bit. No empty site needed.

M1 collapses "is this insert multi-copy in its own genome" to a bit, and throws
away the part that answers ACTIVE. Copy number is a CUMULATIVE record, not a
rate: 30 ISAba1 copies may all be ancient and fixed. It is also family-level, so
it cannot speak to whether THIS insertion is recent.

The event-level, recency-bearing quantity is the focal copy's divergence from
its nearest kin in the SAME genome. An IS copy starts drifting the moment it
lands, so near-identity to its closest relative means the duplication was late.
That needs no empty site and no second genome.

  n_copies              copies at >=--copy-ident identity, >=--copy-cov coverage
  n_distinct_contexts   copies whose +/-200 bp flanks are mutually non-homologous
  min_div_to_kin        NOT MEASURED -- see the failure note below. Do not
                        publish this field.
  [full_length_frac     REMOVED 2026-09-18 -- CIRCULAR. A "copy" is already
                        defined as >=--copy-ident / >=--copy-cov, so a
                        length-outlier is excluded BEFORE it can be counted as a
                        copy. Measured median 1.000 with 0 of 3000 below 0.9: it
                        cannot fire by construction, it is not evidence that
                        decayed copies are absent.]

  CONSEQUENCE that must travel with the catalogue: **n_copies is a count of
  INTACT copies, not of all copies.** Decayed, pseudogenised copies are excluded
  by definition. For an activity estimate that is the right choice, but IS
  annotation tools count every copy, so these numbers read LOW against them.
  Unstated, that difference will be mistaken for disagreement.

**n_distinct_contexts >= 2 is the strongest transposition evidence obtainable
from a single genome**: one element occupying several non-homologous targets
records several independent insertions. n_copies >= 2 with ONE context is most
likely a segmental duplication that carried the element along -- M1 currently
merges these two cases.

PRECONDITION, and it is load-bearing: this works only because the census admits
COMPLETE + CHROMOSOME assemblies. Short-read assemblies collapse identical
repeats, so IS copy number is systematically under-counted at contig level.
These fields become invalid the moment the pool admits draft genomes.

FAILURE NOTE -- min_div_to_kin, 2026-09-18. It did not vary with carrier
frequency in either species (median 0.00000 in every bin). Recorded as NOT
MEASURED, not as measured-zero, and there are TWO independent causes; fixing one
does not settle it:

  A  RESOLUTION. minimap2's NM/blen saturates above ~98% identity, so 92.9%
     (A. baumannii) and 81.8% (MTB) of values are exactly 0.0, max 0.019. Needs
     base-level alignment or mismatch counting from the CIGAR.
  B  THE STATISTIC. A truly recent duplication IS literally identical, so zero
     is the CORRECT answer for it. If so, no resolution rescues a per-bin
     MEDIAN: the age signal lives in the tail. The right statistics would be the
     FRACTION of copies with non-zero divergence, the focal copy's deviation
     from the family consensus, or focal deviation over maximum within-family
     divergence (family age).

  A and B are not separable until A is fixed. Both can be fatal.

AND THE TRUTH AXIS MAY ALSO BE BROKEN. Carrier frequency confounds clonal
expansion with age: a locus at high frequency may be old, or may be recent in a
lineage that then expanded. So the flat result implicates the instrument, the
statistic, OR the ruler. Practice 7 (invariance is a warning) points in more
than one direction here.

    107_within_genome_copies.py --recon <dirs> --events <tsv> \\
        --genome-dir <dir> --species x --out x
"""
import argparse, collections, csv, glob, os, random, statistics, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


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


def load_genome(p):
    seqs, cur = {}, None
    for line in open(p):
        if line[0] == ">":
            cur = line[1:].split()[0]; seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--copy-ident", type=float, default=0.98)
    ap.add_argument("--copy-cov", type=float, default=0.90)
    ap.add_argument("--flank", type=int, default=200)
    ap.add_argument("--flank-ident", type=float, default=0.90)
    ap.add_argument("--sample", type=int, default=3000)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    recon = []
    for d in args.recon:
        recon.extend(sorted(glob.glob(d)) if any(c in d for c in "*?[") else [d])
    loci, aseq, carr, alln = {}, {}, {}, {}
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
                g = [x for x in r.get("genomes", "").split(",") if x]
                carr.setdefault(r["locus_id"], {})[r["allele_id"]] = g
                alln[r["locus_id"]] = alln.get(r["locus_id"], 0) + len(g)
    ev = list(csv.DictReader(open(args.events), delimiter="\t"))
    log("%d events, %d loci loaded" % (len(ev), len(loci)))

    random.seed(3)
    sel = random.sample(ev, min(args.sample, len(ev)))
    bygen = collections.defaultdict(list)
    for r in sel:
        lid = r["locus_id"]
        row = loci.get(lid)
        if not row:
            continue
        g = (carr.get(lid, {}).get(row["longest_allele"]) or [None])[0]
        if g:
            bygen[g].append(lid)
    log("%d events map to %d carrier genomes" % (sum(len(v) for v in bygen.values()), len(bygen)))

    tmp = tempfile.mkdtemp(prefix="wgc.")
    out = {}
    for g, lids in bygen.items():
        fa = os.path.join(args.genome_dir, g + ".fna")
        if not os.path.exists(fa):
            continue
        q = os.path.join(tmp, "q.fa")
        with open(q, "w") as fh:
            for lid in lids:
                row = loci[lid]
                ls = aseq.get(lid, {}).get(row["longest_allele"], "")
                o, il = int(row["lcp_bp"]), int(row["inserted_len"])
                if ls:
                    fh.write(">%s\n%s\n" % (lid, ls[o:o + il]))
        paf = os.path.join(tmp, "p.paf")
        with open(paf, "w") as fh:
            subprocess.run(["minimap2", "-c", "-x", "asm10", "-N", "80", "-p", "0.05",
                            "-t", str(args.threads), fa, q], stdout=fh,
                           stderr=subprocess.DEVNULL, check=False)
        hits = collections.defaultdict(list)
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qn, ql = f[0], int(f[1])
            nm, bl = int(f[9]), int(f[10])
            if bl == 0 or ql == 0:
                continue
            idt, cov = nm / bl, bl / ql
            if idt >= args.copy_ident and cov >= args.copy_cov:
                hits[qn].append((f[5], int(f[7]), int(f[8]), idt, bl))
        gen = load_genome(fa)
        for lid in lids:
            h = hits.get(lid, [])
            if not h:
                out[lid] = (0, 0, float("nan"), float("nan"))
                continue
            # distinct contexts: single-linkage on flank similarity
            flanks = []
            for c, s, e, idt, bl in h:
                seq = gen.get(c, "")
                flanks.append((seq[max(0, s - args.flank):s],
                               seq[e:e + args.flank]))
            ctx, used = 0, [False] * len(h)
            for i in range(len(h)):
                if used[i]:
                    continue
                ctx += 1
                used[i] = True
                for j in range(i + 1, len(h)):
                    if used[j]:
                        continue
                    same = 0
                    for a_, b_ in ((flanks[i][0], flanks[j][0]),
                                   (flanks[i][1], flanks[j][1]),
                                   (flanks[i][0], rc(flanks[j][1])),
                                   (flanks[i][1], rc(flanks[j][0]))):
                        if a_ and b_ and len(a_) > 50 and len(b_) > 50:
                            m = sum(1 for x, y in zip(a_, b_) if x == y)
                            if m / min(len(a_), len(b_)) >= args.flank_ident:
                                same = 1
                                break
                    if same:
                        used[j] = True
            idents = sorted((x[3] for x in h), reverse=True)
            mind = 1.0 - idents[1] if len(idents) > 1 else float("nan")
            out[lid] = (len(h), ctx, mind, float("nan"))

    cols = ["event_id", "locus_id", "species", "inserted_len", "M1_within_genome",
            "n_copies", "n_distinct_contexts", "min_div_to_kin", "full_length_frac",
            "carrier_frac"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    rows = []
    for r in sel:
        lid = r["locus_id"]
        if lid not in out:
            continue
        n, c, d_, f_ = out[lid]
        row = loci.get(lid, {})
        nc = len(carr.get(lid, {}).get(row.get("longest_allele", ""), []))
        tot = alln.get(lid, 0)
        rows.append({"event_id": r["event_id"], "locus_id": lid,
                     "species": args.species, "inserted_len": r["inserted_len"],
                     "M1_within_genome": r.get("M1_within_genome", "."),
                     "n_copies": n, "n_distinct_contexts": c,
                     "min_div_to_kin": "%.5f" % d_ if d_ == d_ else "NA",
                     "full_length_frac": "%.3f" % f_ if f_ == f_ else "NA",
                     "carrier_frac": "%.3f" % (nc / tot) if tot else "NA"})
    with open(args.out + "_wgcopies.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    log("=" * 76)
    log("WITHIN-GENOME COPY STRUCTURE -- %s   (n=%d)" % (args.species, len(rows)))
    log("")
    h = collections.Counter()
    for r in rows:
        n, c = r["n_copies"], r["n_distinct_contexts"]
        h["0-1 copy" if n <= 1 else
          ">=2 copies, 1 context (segmental dup?)" if c <= 1 else
          ">=2 copies, >=2 CONTEXTS (transposition)"] += 1
    for k in ("0-1 copy", ">=2 copies, 1 context (segmental dup?)",
              ">=2 copies, >=2 CONTEXTS (transposition)"):
        log("  %-42s %6d  %5.1f%%" % (k, h[k], 100.0 * h[k] / len(rows) if rows else 0))
    log("")
    log("  M1 merges the middle and bottom rows; they are different events.")
    log("")
    log("  FALSIFICATION -- does min_div_to_kin fall with carrier frequency?")
    log("  (a locus in 1 of 12 genomes is younger than one in 11 of 12)")
    log("  %-14s %8s %14s" % ("carrier_frac", "n", "med div_to_kin"))
    log("  " + "-" * 38)
    for lo, hi, lab in ((0, .2, "<0.2"), (.2, .4, "0.2-0.4"), (.4, .6, "0.4-0.6"),
                        (.6, .8, "0.6-0.8"), (.8, 1.01, ">=0.8")):
        v = [float(r["min_div_to_kin"]) for r in rows
             if r["min_div_to_kin"] != "NA" and r["carrier_frac"] != "NA"
             and lo <= float(r["carrier_frac"]) < hi]
        if v:
            log("  %-14s %8d %14.5f" % (lab, len(v), statistics.median(v)))
    log("")
    log("  If this does not RISE with carrier frequency, min_div_to_kin is not")
    log("  measuring recency and the field must not be read as one.")
    log("=" * 76)
    return 0


if __name__ == "__main__":
    sys.exit(main())
