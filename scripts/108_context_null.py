#!/usr/bin/env python3
"""Composition-matched null for `n_distinct_contexts >= 2`. Run BEFORE publishing it.

`>=2 distinct contexts` is read as "one element at several non-homologous
targets in one genome = multiple independent insertions". A class of natural
false positives passes it perfectly WITHOUT ever transposing: rRNA operons,
multi-copy gene families, phage genes -- the same sequence in several
non-homologous surroundings, arrived by other means.

So the statistic needs a baseline (practice 2: score against a measured,
composition-matched null, never a bare threshold). The null here is a random
genomic interval from the SAME carrier genome, matched on LENGTH and GC, not a
bubble, put through the identical counter.

  baseline ~40%  -> M. tuberculosis's 42.3% says nothing; A. baumannii's 88.9%
                    still stands
  baseline low   -> both species' figures stand

KNOWN AND UNFIXED, a second-order boundary: two independent insertions into the
same hotspot share flanking sequence, so the 200 bp flank-homology test merges
them into ONE context. That is a systematic UNDER-count, and it falls hardest on
exactly the loci a target-site study cares about. Noted, not corrected.

    108_context_null.py --recon <dirs> --events <tsv> --genome-dir <dir> \\
        --species x --out x
"""
import argparse, collections, csv, glob, os, random, statistics, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def gc(s):
    n = sum(1 for c in s if c in "GC")
    return n / len(s) if s else 0.0


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


def count_contexts(q_records, fa, gen, args):
    """Identical counter for real inserts and null intervals."""
    tmp = tempfile.mkdtemp(prefix="ctxnull.")
    q = os.path.join(tmp, "q.fa")
    with open(q, "w") as fh:
        for name, seq in q_records:
            fh.write(">%s\n%s\n" % (name, seq))
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
        if nm / bl >= args.copy_ident and bl / ql >= args.copy_cov:
            hits[qn].append((f[5], int(f[7]), int(f[8])))
    res = {}
    for name, _ in q_records:
        h = hits.get(name, [])
        if not h:
            res[name] = (0, 0)
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
                            used[j] = True
                            break
        res[name] = (len(h), ctx)
    return res


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
    ap.add_argument("--gc-tol", type=float, default=0.03)
    ap.add_argument("--sample", type=int, default=800)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    recon = []
    for d in args.recon:
        recon.extend(sorted(glob.glob(d)) if any(c in d for c in "*?[") else [d])
    loci, aseq, carr = {}, {}, {}
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
                carr.setdefault(r["locus_id"], {})[r["allele_id"]] = \
                    [x for x in r.get("genomes", "").split(",") if x]
    ev = list(csv.DictReader(open(args.events), delimiter="\t"))
    random.seed(21)
    sel = random.sample(ev, min(args.sample, len(ev)))
    bygen = collections.defaultdict(list)
    for r in sel:
        lid = r["locus_id"]; row = loci.get(lid)
        if not row:
            continue
        g = (carr.get(lid, {}).get(row["longest_allele"]) or [None])[0]
        if g:
            bygen[g].append(lid)
    log("%d events across %d carrier genomes" % (sum(len(v) for v in bygen.values()), len(bygen)))

    R, N = collections.Counter(), collections.Counter()
    nreal = nnull = 0
    for g, lids in bygen.items():
        fa = os.path.join(args.genome_dir, g + ".fna")
        if not os.path.exists(fa):
            continue
        gen = load_genome(fa)
        big = [c for c in gen if len(gen[c]) > 100000]
        if not big:
            continue
        real, null = [], []
        for lid in lids:
            row = loci[lid]
            ls = aseq.get(lid, {}).get(row["longest_allele"], "")
            o, il = int(row["lcp_bp"]), int(row["inserted_len"])
            if not ls:
                continue
            ins = ls[o:o + il]
            real.append((lid, ins))
            # length- and GC-matched random interval from the same genome
            tgt = gc(ins)
            for _ in range(40):
                c = random.choice(big)
                st = random.randrange(0, len(gen[c]) - il)
                cand = gen[c][st:st + il]
                if cand.count("N") > 0.01 * il:
                    continue
                if abs(gc(cand) - tgt) <= args.gc_tol:
                    null.append(("NULL_" + lid, cand))
                    break
        if real:
            rr = count_contexts(real, fa, gen, args)
            for k, (n, c) in rr.items():
                nreal += 1
                R[">=2ctx" if c >= 2 else "1ctx" if n >= 2 else "0-1copy"] += 1
        if null:
            nn = count_contexts(null, fa, gen, args)
            for k, (n, c) in nn.items():
                nnull += 1
                N[">=2ctx" if c >= 2 else "1ctx" if n >= 2 else "0-1copy"] += 1

    log("=" * 74)
    log("CONTEXT NULL -- %s" % args.species)
    log("  null = random interval, SAME genome, matched length and GC (+/-%.2f)"
        % args.gc_tol)
    log("")
    log("  %-28s %10s %10s" % ("", "REAL insert", "NULL interval"))
    log("  " + "-" * 50)
    for k in ("0-1copy", "1ctx", ">=2ctx"):
        lab = {"0-1copy": "0-1 copy", "1ctx": ">=2 copies, 1 context",
               ">=2ctx": ">=2 copies, >=2 CONTEXTS"}[k]
        log("  %-28s %9.1f%% %9.1f%%" % (lab,
            100.0 * R[k] / nreal if nreal else 0,
            100.0 * N[k] / nnull if nnull else 0))
    log("")
    log("  n real=%d  n null=%d" % (nreal, nnull))
    rr = 100.0 * R[">=2ctx"] / nreal if nreal else 0
    nn_ = 100.0 * N[">=2ctx"] / nnull if nnull else 0
    log("  >=2 contexts: real %.1f%%  vs  null %.1f%%   ratio %.2fx"
        % (rr, nn_, rr / nn_ if nn_ else float("inf")))
    if nn_ >= 0.5 * rr:
        log("  *** THE BASELINE EXPLAINS MOST OF IT. Do not read >=2 contexts as")
        log("  transposition evidence for this species. ***")
    log("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
