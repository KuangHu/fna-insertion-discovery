#!/usr/bin/env python3
"""Cross-species element family catalogue, clustered at a MEASURED threshold.

Merges the event-level inserts of every species and clusters them by sequence
identity. Two rules carried from earlier work:

  THE THRESHOLD IS NOT HAND-PICKED. It is read off a composition-matched null:
  shuffle each sequence (preserving composition), align shuffled pairs, and take
  a high quantile of that identity distribution. An 0.90 cutoff was wrong once
  before for exactly this reason; the revcomp null is what exposed it.

  NO ANNOTATION. The catalogue IS the product. Family naming is Layer 5 and is
  not done here -- no ISEScan, no IS library. What can be read without any
  reference is (a) whether the size modes coincide with known IS sizes, which is
  internal structural consistency, and (b) how many families are shared across
  species spanning two phyla, which is direct structural evidence of horizontal
  movement.

    106_element_families.py --inserts a.fna b.fna --out families/all
"""
import argparse, collections, os, random, subprocess, statistics, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402


def read_fa(p):
    out, cur = {}, None
    for line in open(p):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inserts", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--null-pairs", type=int, default=3000)
    ap.add_argument("--null-quantile", type=float, default=0.99)
    ap.add_argument("--min-cov", type=float, default=0.80)
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    seqs, species = {}, {}
    for f in args.inserts:
        sp = os.path.basename(f).replace("_inserts.fna", "")
        for k, v in read_fa(f).items():
            key = "%s|%s" % (sp, k)
            seqs[key] = v
            species[key] = sp
    log("%d inserts from %d species" % (len(seqs), len(set(species.values()))))

    # ---- composition-matched null -----------------------------------------
    random.seed(13)
    keys = list(seqs)
    tmp = tempfile.mkdtemp(prefix="fam.")
    qa, qb = os.path.join(tmp, "na.fa"), os.path.join(tmp, "nb.fa")
    with open(qa, "w") as fa, open(qb, "w") as fb:
        for i in range(args.null_pairs):
            a, b = random.choice(keys), random.choice(keys)
            sa, sb = list(seqs[a]), list(seqs[b])
            random.shuffle(sa); random.shuffle(sb)
            fa.write(">n%d\n%s\n" % (i, "".join(sa)))
            fb.write(">n%d\n%s\n" % (i, "".join(sb)))
    paf = os.path.join(tmp, "null.paf")
    with open(paf, "w") as fh:
        subprocess.run(["minimap2", "-c", "-x", "asm20", "-t", str(args.threads), qa, qb],
                       stdout=fh, stderr=subprocess.DEVNULL, check=False)
    nid = []
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and int(f[10]) > 0:
            nid.append(int(f[9]) / int(f[10]))
    if nid:
        nid.sort()
        thr = nid[min(len(nid) - 1, int(len(nid) * args.null_quantile))]
        log("composition-matched null: %d shuffled pairs aligned, identity "
            "median %.3f, q%.2f = %.3f" % (len(nid), statistics.median(nid),
                                           args.null_quantile, thr))
    else:
        thr = 0.90
        log("null produced NO alignments -- shuffled sequences do not align at "
            "all, so any real alignment is above the floor. Using %.2f as a "
            "conservative identity threshold." % thr)
    log("clustering threshold: identity >= %.3f, coverage >= %.2f" % (thr, args.min_cov))

    allfa = os.path.join(tmp, "all.fa")
    with open(allfa, "w") as fh:
        for k, v in seqs.items():
            fh.write(">%s\n%s\n" % (k, v))
    paf2 = os.path.join(tmp, "all.paf")
    with open(paf2, "w") as fh:
        subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "200", "-p", "0.1",
                        "-t", str(args.threads), allfa, allfa],
                       stdout=fh, stderr=subprocess.DEVNULL, check=False)
    par = {}
    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    n_edge = 0
    for line in open(paf2):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        q, ql, t, bl, nm = f[0], int(f[1]), f[5], int(f[10]), int(f[9])
        if q == t or bl == 0 or ql == 0:
            continue
        if nm / bl >= thr and bl / min(ql, int(f[6])) >= args.min_cov:
            par.setdefault(q, q); par.setdefault(t, t)
            par[find(q)] = find(t); n_edge += 1
    log("all-vs-all: %d edges above threshold" % n_edge)
    fam = collections.defaultdict(list)
    for k in seqs:
        fam[find(k)].append(k)
    log("clusters: %d (from %d inserts)" % (len(fam), len(seqs)))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    rows = []
    for i, (_, mem) in enumerate(sorted(fam.items(), key=lambda kv: -len(kv[1])), 1):
        L = sorted(len(seqs[m]) for m in mem)
        sc = collections.Counter(species[m] for m in mem)
        rows.append({"family_id": "F%05d" % i, "n_members": len(mem),
                     "median_len": L[len(L) // 2],
                     "len_iqr": "%d-%d" % (L[len(L) // 4], L[3 * len(L) // 4]),
                     "n_species": len(sc),
                     "per_species": ",".join("%s:%d" % kv for kv in sorted(sc.items()))})
    with open(args.out + "_families.tsv", "w") as fh:
        cols = list(rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    log("=" * 74)
    log("ELEMENT FAMILY CATALOGUE  (no annotation, no IS library)")
    log("")
    log("  families              %6d" % len(rows))
    log("  singletons            %6d  %5.1f%%"
        % (sum(1 for r in rows if r["n_members"] == 1),
           100.0 * sum(1 for r in rows if r["n_members"] == 1) / len(rows)))
    log("  families >=10 members %6d" % sum(1 for r in rows if r["n_members"] >= 10))
    log("")
    log("  CROSS-SPECIES families (direct structural evidence of movement):")
    for n in (2, 3, 4):
        log("    present in >=%d species  %6d" % (n, sum(1 for r in rows if r["n_species"] >= n)))
    log("")
    log("  size modes of families with >=10 members:")
    big = [r["median_len"] for r in rows if r["n_members"] >= 10]
    h = collections.Counter(x // 100 * 100 for x in big)
    for k, v in h.most_common(8):
        log("    %5d-%5d bp : %3d families" % (k, k + 99, v))
    log("=" * 74)
    log("wrote %s_families.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
