#!/usr/bin/env python3
"""Why did Arm A miss these IS? A taxonomy, NOT a parameter sweep.

Stage 60 says Arm A misses 502 of 1086 multi-copy ISEScan calls. That number on
its own cannot tell you whether to change the clustering, lower the identity
threshold, or accept the regime -- so this script classifies every miss by the
point in Arm A at which it was lost, and stops there. No thresholds are tuned
here and none should be until this table is read.

The decisive question is whether minimap2's self-alignment SAW the repeat at
all. That splits every miss into two fundamentally different populations:

  * self-alignment never produced an interval at the locus
      -> the repeat is invisible to Arm A's front end. Lowering a downstream
         threshold cannot recover it. Causes: copies too diverged to align,
         copies collapsed by the assembler, or a single-copy IS in this genome
         that ISEScan called from its pHMM rather than from copy structure.
  * self-alignment DID produce an interval, and something later dropped it
      -> recoverable. The sub-class names which gate did the dropping.

Classes, applied in this order (first match wins):

  OUT_OF_SIZE_WINDOW      IS length outside [min_len, max_len]; by design
  CONTIG_EDGE             IS within --flank of a contig end; QC-excluded by
                          design, cannot be repaired by algorithm
  EXTENT_MISMATCH_*       Arm A HAS an interval here but reciprocal overlap
                          < 0.5, so stage 60 scored it a miss. This is a
                          BOUNDARY disagreement, not a discovery failure, and
                          must not be pooled with the real misses.
                            _ARMA_LONGER   merged/nested repeat family
                            _ARMA_SHORTER  truncated call
  SELF_ALN_ABSENT         no self-alignment interval covers the locus
  BELOW_MIN_COPIES        self-aln found it, but < --min-copies distinct loci
  BELOW_MIN_IDENTITY      self-aln found it, best copy identity < --min-identity
  TANDEM_ONLY             all copies within --min-locus-sep: a tandem
                          duplication, which Arm A excludes on purpose
  FILTERED_DOWNSTREAM     passed every gate above and was still lost, so the
                          loss is in boundary voting / sharpness / flank
                          homology -- the genuinely interesting bucket

    62_armA_miss_taxonomy.py --per-is bench60_per_is.tsv \\
        --paf-glob 'armA_scale/*_self.paf' --fna-dir <dir> --out miss_tax
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log, overlap                               # noqa: E402


def load_paf(path):
    """Self-alignment intervals, keyed by query contig. Diagonal hits dropped."""
    by_contig = defaultdict(list)
    clen = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            t, tl, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
            nmatch, alen = int(f[9]), max(1, int(f[10]))
            clen[q], clen[t] = ql, tl
            if q == t and qs == ts and qe == te:
                continue                       # trivial self hit
            by_contig[q].append({"qs": qs, "qe": qe, "t": t, "ts": ts,
                                 "te": te, "ident": nmatch / alen})
    return by_contig, clen


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-is", required=True, help="bench60_per_is.tsv")
    ap.add_argument("--paf-glob", required=True, help="'*_self.paf'")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=300)
    ap.add_argument("--max-len", type=int, default=5000)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--min-copies", type=int, default=3)
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-locus-sep", type=int, default=10000)
    ap.add_argument("--cover-frac", type=float, default=0.5,
                    help="fraction of the IS an interval must cover to count")
    args = ap.parse_args()

    pafs = {}
    for p in sorted(glob.glob(args.paf_glob)):
        s = os.path.basename(p)[: -len("_self.paf")]
        pafs[s] = p
    log("%d self-alignment PAFs available" % len(pafs))

    cache = {}

    def get(sample):
        if sample not in cache:
            cache[sample] = (load_paf(pafs[sample]) if sample in pafs
                             else ({}, {}))
        return cache[sample]

    rows, counts = [], defaultdict(int)
    by_family = defaultdict(lambda: defaultdict(int))
    with open(args.per_is) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["class"] != "multi" or r["found_by_armA"] != "0":
                continue
            s, c = r["sample"], r["contig"]
            a, b = int(r["is_start"]), int(r["is_end"])
            L = b - a
            ov = float(r["reciprocal_overlap"] or 0)
            by_contig, clen = get(s)

            hits = [h for h in by_contig.get(c, ())
                    if overlap(a, b, h["qs"], h["qe"]) >= args.cover_frac * L]
            n_loci, best_id, seps = 0, 0.0, []
            if hits:
                loci = {(h["t"], h["ts"] // 1000) for h in hits}
                n_loci = len(loci) + 1              # + the query locus itself
                best_id = max(h["ident"] for h in hits)
                seps = [abs(h["ts"] - a) if h["t"] == c else args.min_locus_sep + 1
                        for h in hits]

            cl = clen.get(c, 0)
            edge = min(a, cl - b) if cl else args.flank + 1

            if L < args.min_len or L > args.max_len:
                k = "OUT_OF_SIZE_WINDOW"
            elif cl and edge < args.flank:
                k = "CONTIG_EDGE"
            elif ov > 0:
                k = ("EXTENT_MISMATCH_ARMA_LONGER" if ov < 0.5 and r["armA_family"] != "."
                     else "EXTENT_MISMATCH")
            elif not hits:
                k = "SELF_ALN_ABSENT"
            elif n_loci < args.min_copies:
                k = "BELOW_MIN_COPIES"
            elif best_id < args.min_identity:
                k = "BELOW_MIN_IDENTITY"
            elif seps and max(seps) < args.min_locus_sep:
                k = "TANDEM_ONLY"
            else:
                k = "FILTERED_DOWNSTREAM"

            counts[k] += 1
            by_family[r["is_family"]][k] += 1
            rows.append([s, c, a, b, L, r["is_family"], r["is_copy_number"],
                         "%.3f" % ov, n_loci, "%.4f" % best_id, edge, k])

    with open(args.out + "_per_miss.tsv", "w") as fh:
        fh.write("sample\tcontig\tis_start\tis_end\tis_len\tis_family\t"
                 "is_copy_number\treciprocal_overlap\tn_selfaln_loci\t"
                 "best_selfaln_identity\tdist_to_contig_end\tmiss_class\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    tot = sum(counts.values())
    with open(args.out + "_summary.tsv", "w") as fh:
        fh.write("miss_class\tn\tpct\n")
        for k, v in sorted(counts.items(), key=lambda x: -x[1]):
            fh.write("%s\t%d\t%.1f\n" % (k, v, 100.0 * v / tot if tot else 0))

    with open(args.out + "_by_family.tsv", "w") as fh:
        classes = sorted(counts)
        fh.write("is_family\ttotal\t" + "\t".join(classes) + "\n")
        for fam in sorted(by_family, key=lambda k: -sum(by_family[k].values())):
            d = by_family[fam]
            fh.write("%s\t%d\t%s\n" % (fam, sum(d.values()),
                                       "\t".join(str(d.get(c, 0)) for c in classes)))

    log("=" * 72)
    log("ARM A MISS TAXONOMY  (n=%d multi-copy ISEScan calls scored as missed)" % tot)
    log("")
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        log("  %-32s %5d  %5.1f%%" % (k, v, 100.0 * v / tot if tot else 0))
    log("")
    recoverable = sum(v for k, v in counts.items()
                      if k in ("BELOW_MIN_COPIES", "BELOW_MIN_IDENTITY",
                               "FILTERED_DOWNSTREAM"))
    bydesign = sum(v for k, v in counts.items()
                   if k in ("OUT_OF_SIZE_WINDOW", "CONTIG_EDGE", "TANDEM_ONLY"))
    scoring = sum(v for k, v in counts.items() if k.startswith("EXTENT_MISMATCH"))
    blind = counts.get("SELF_ALN_ABSENT", 0)
    log("  recoverable by threshold/logic change : %5d  %5.1f%%"
        % (recoverable, 100.0 * recoverable / tot if tot else 0))
    log("  excluded by design (not a bug)        : %5d  %5.1f%%"
        % (bydesign, 100.0 * bydesign / tot if tot else 0))
    log("  boundary disagreement, not a miss     : %5d  %5.1f%%"
        % (scoring, 100.0 * scoring / tot if tot else 0))
    log("  invisible to the self-alignment       : %5d  %5.1f%%"
        % (blind, 100.0 * blind / tot if tot else 0))
    log("=" * 72)
    log("wrote %s_{per_miss,summary,by_family}.tsv" % args.out)


if __name__ == "__main__":
    main()
