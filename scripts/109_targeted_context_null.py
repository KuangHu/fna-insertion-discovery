#!/usr/bin/env python3
"""TARGETED null for `>=2 distinct contexts`: multi-copy sequences that are NOT insertions.

The random-interval null (108_) answered the wrong question. Its 172x / 291x
ratios come almost entirely from the fact that a random piece of a genome is not
repetitive at all -- 97.7% and 99.6% of random intervals have 0-1 copies. The
worry was never that. It was: **do sequences that are ALREADY multi-copy --
rRNA operons, multi-copy gene families, phage genes -- occupy several
non-homologous contexts without ever having transposed?**

Conditioning the random null on n_copies >= 2 answers that, but the denominator
collapses: A. baumannii 4/18, M. tuberculosis **1/3**. Unusable, and MTB's is
nothing at all.

So condition by CONSTRUCTION instead. Arm A's families are multi-copy by its own
criterion (>=3 copies at distinct loci), found by self-alignment with no
insertion model. Remove the ones that correspond to a called insertion, and what
remains is exactly the target population: multi-copy, not called as an
insertion. Run the identical counter over it.

  real insertions   >=2 contexts at rate R
  non-insertion
  multi-copy seqs   >=2 contexts at rate N     <- the number that matters

WHAT THIS STILL DOES NOT MATCH, recorded rather than fixed: bubble sequences
come by construction from regions that DIFFER between genomes, i.e. the
accessory genome, which is more repeat-rich than core. Arm A families are drawn
from the whole genome. So the comparison remains somewhat generous. The tighter
null is other BUBBLES that are not candidate insertions (SNP/indel-type,
anchor-less); left for re-checking.

    109_targeted_context_null.py --arma <dir> --inserts <fna> \\
        --genome-dir <dir> --species x --out x
"""
import argparse, collections, csv, glob, os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(p):
    out, cur = {}, None
    if not os.path.exists(p):
        return out
    for line in open(p):
        if line.startswith(">"):
            cur = line[1:].split()[0]; out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arma", required=True)
    ap.add_argument("--inserts", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--copy-ident", type=float, default=0.98)
    ap.add_argument("--copy-cov", type=float, default=0.90)
    ap.add_argument("--flank", type=int, default=200)
    ap.add_argument("--flank-ident", type=float, default=0.90)
    ap.add_argument("--is-insert-ident", type=float, default=0.90)
    ap.add_argument("--sample", type=int, default=600)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="tgtnull.")
    # 1. every Arm A family element, grouped by its genome
    bygen = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(args.arma, "*_elements.fna")):
        g = os.path.basename(f)[:-len("_elements.fna")]
        for k, v in read_fa(f).items():
            bygen[g][k] = v
    n_el = sum(len(v) for v in bygen.values())
    log("%d Arm A family elements across %d genomes" % (n_el, len(bygen)))
    if n_el == 0:
        log("FATAL: no Arm A elements found")
        return 2

    # 2. drop the ones that ARE a called insertion
    allel = os.path.join(tmp, "el.fa")
    with open(allel, "w") as fh:
        for g, d in bygen.items():
            for k, v in d.items():
                fh.write(">%s\n%s\n" % (k, v))
    paf = os.path.join(tmp, "ins.paf")
    with open(paf, "w") as fh:
        subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1",
                        "-t", str(args.threads), allel, args.inserts],
                       stdout=fh, stderr=subprocess.DEVNULL, check=False)
    is_ins = set()
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        nm, bl, tl = int(f[9]), int(f[10]), int(f[6])
        if bl and nm / bl >= args.is_insert_ident and bl / tl >= 0.5:
            is_ins.add(f[5])
    log("of those, %d match a called insertion (>=%.2f id, >=50%% of the element) "
        "-> excluded" % (len(is_ins), args.is_insert_ident))
    pool = [(g, k, v) for g, d in bygen.items() for k, v in d.items() if k not in is_ins]
    log("TARGETED NULL POOL: %d multi-copy, non-insertion sequences" % len(pool))
    if len(pool) < 30:
        log("FATAL: pool of %d is too small to estimate a rate" % len(pool))
        return 2
    random.seed(31)
    sel = random.sample(pool, min(args.sample, len(pool)))
    bysel = collections.defaultdict(list)
    for g, k, v in sel:
        bysel[g].append((k, v))

    N = collections.Counter(); tot = 0
    for g, recs in bysel.items():
        fa = os.path.join(args.genome_dir, g + ".fna")
        if not os.path.exists(fa):
            continue
        gen = read_fa(fa)
        q = os.path.join(tmp, "q.fa")
        with open(q, "w") as fh:
            for k, v in recs:
                fh.write(">%s\n%s\n" % (k, v))
        p2 = os.path.join(tmp, "p.paf")
        with open(p2, "w") as fh:
            subprocess.run(["minimap2", "-c", "-x", "asm10", "-N", "80", "-p", "0.05",
                            "-t", str(args.threads), fa, q], stdout=fh,
                           stderr=subprocess.DEVNULL, check=False)
        hits = collections.defaultdict(list)
        for line in open(p2):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qn, ql = f[0], int(f[1]); nm, bl = int(f[9]), int(f[10])
            if bl and ql and nm / bl >= args.copy_ident and bl / ql >= args.copy_cov:
                hits[qn].append((f[5], int(f[7]), int(f[8])))
        for k, _ in recs:
            h = hits.get(k, [])
            tot += 1
            if len(h) < 2:
                N["0-1copy"] += 1
                continue
            fl = []
            for c, s, e in h:
                sq = gen.get(c, "")
                fl.append((sq[max(0, s - args.flank):s], sq[e:e + args.flank]))
            ctx, used = 0, [False] * len(h)
            for i in range(len(h)):
                if used[i]:
                    continue
                ctx += 1; used[i] = True
                for j in range(i + 1, len(h)):
                    if used[j]:
                        continue
                    for a_, b_ in ((fl[i][0], fl[j][0]), (fl[i][1], fl[j][1]),
                                   (fl[i][0], rc(fl[j][1])), (fl[i][1], rc(fl[j][0]))):
                        if a_ and b_ and len(a_) > 50 and len(b_) > 50:
                            m = sum(1 for x, y in zip(a_, b_) if x == y)
                            if m / min(len(a_), len(b_)) >= args.flank_ident:
                                used[j] = True; break
            N[">=2ctx" if ctx >= 2 else "1ctx"] += 1

    log("=" * 74)
    log("TARGETED CONTEXT NULL -- %s" % args.species)
    log("  population: Arm A multi-copy families that are NOT called insertions")
    log("")
    for k in ("0-1copy", "1ctx", ">=2ctx"):
        lab = {"0-1copy": "0-1 copy (re-mapped)", "1ctx": ">=2 copies, 1 context",
               ">=2ctx": ">=2 copies, >=2 CONTEXTS"}[k]
        log("  %-30s %6d  %5.1f%%" % (lab, N[k], 100.0 * N[k] / tot if tot else 0))
    multi = N["1ctx"] + N[">=2ctx"]
    log("")
    log("  n scored = %d ;  multi-copy among them = %d" % (tot, multi))
    if multi >= 30:
        log("  CONDITIONAL on multi-copy: %d/%d = %.1f%% span >=2 contexts"
            % (N[">=2ctx"], multi, 100.0 * N[">=2ctx"] / multi))
        log("  (real insertions: A. baumannii 98.8%%, M. tuberculosis 95.2%%)")
    else:
        log("  CONDITIONAL NOT REPORTED: only %d multi-copy sequences, too few" % multi)
    log("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
