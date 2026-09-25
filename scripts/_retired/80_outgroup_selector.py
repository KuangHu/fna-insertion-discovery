#!/usr/bin/env python3
"""Layer 3a -- choose outgroup genomes for a cluster. Selection only, no polarity.

Polarity needs genomes OUTSIDE the near-identical cluster, but "outside" is a
band, not a direction: too close and the outgroup shares the event under test,
too far and Layer 2's 500 bp anchors stop placing and the locus cannot be
compared at all. So an outgroup must clear four independent conditions:

  1. not a member of the ingroup cluster;
  2. whole-genome relatedness inside a band (--min-ani/--max-ani), high enough
     that the locus is still homologous;
  3. aligned fraction high enough that anchors have somewhere to land
     (--min-af) -- ANI alone is not sufficient, two genomes can be 99%
     identical over the 30% of themselves that aligns;
  4. INDEPENDENCE from the other chosen outgroups. Three near-identical
     outgroups are one observation wearing three hats, and would inflate any
     later agreement fraction. Candidates within --max-outgroup-ani of an
     already-selected one are skipped.

Nothing here decides which allele is ancestral, or even looks at a locus. This
script answers only "which genomes are entitled to vote later".

    80_outgroup_selector.py --ingroup cl0000.manifest --pool-list all_fna.txt \\
        --outdir outgroup/cl0000 --min-ani 95 --max-ani 99 --min-af 80 -n 5
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                        # noqa: E402


def sample_name(path):
    b = os.path.basename(path)
    for suf in (".gz", ".fna", ".fa", ".fasta"):
        if b.endswith(suf):
            b = b[: -len(suf)]
    return b


def skani_dist(ql, rl, out, threads, skani, min_af):
    cmd = [skani, "dist", "--ql", ql, "--rl", rl, "-t", str(threads),
           "-o", out, "--min-af", str(min_af)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    return out


def read_dist(path):
    """(query_sample, ref_sample) -> (ani, af_ref, af_query)"""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, t = sample_name(r["Query_file"]), sample_name(r["Ref_file"])
            try:
                out[(q, t)] = (float(r["ANI"]),
                               float(r["Align_fraction_ref"]),
                               float(r["Align_fraction_query"]))
            except (ValueError, KeyError):
                continue
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ingroup", required=True, help="cluster manifest")
    ap.add_argument("--pool-list", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-ani", type=float, default=95.0)
    ap.add_argument("--max-ani", type=float, default=99.0,
                    help="above this the genome is effectively ingroup and may "
                         "carry the very event whose polarity is in question")
    ap.add_argument("--min-af", type=float, default=80.0,
                    help="aligned fraction; the real guard on anchor placement")
    ap.add_argument("-n", "--n-outgroups", type=int, default=5)
    ap.add_argument("--n-candidates", type=int, default=300,
                    help="how many ranked candidates enter the independence "
                         "pass. Must be MUCH larger than -n: ranking by aligned "
                         "fraction concentrates on one clonal lineage, so a "
                         "short list can dedup down to a single outgroup. "
                         "Measured: with 40, cl0000 yielded 1 outgroup from "
                         "5142 in-band candidates.")
    ap.add_argument("--max-outgroup-ani", type=float, default=99.0,
                    help="selected outgroups must be less similar than this to "
                         "each other, so they are independent observations")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skani", default="skani")
    ap.add_argument("--reuse", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ing = [l.strip() for l in open(args.ingroup) if l.strip()]
    ing_names = {sample_name(p) for p in ing}
    log("ingroup: %d genomes" % len(ing))

    d1 = os.path.join(args.outdir, "ingroup_vs_pool.tsv")
    if not (args.reuse and os.path.exists(d1)):
        log("skani: ingroup vs pool ...")
        skani_dist(args.ingroup, args.pool_list, d1, args.threads, args.skani,
                   min(50.0, args.min_af))
    dist = read_dist(d1)

    # aggregate each candidate over all ingroup members; use the MINIMUM ANI and
    # MINIMUM aligned fraction, so a candidate must be comparable to every
    # ingroup genome, not merely to the most permissive one
    agg = defaultdict(list)
    for (q, t), (ani, afr, afq) in dist.items():
        if t in ing_names:
            continue
        agg[t].append((ani, min(afr, afq)))
    cands = []
    for t, vals in agg.items():
        if len(vals) < len(ing_names):
            continue                       # not comparable to every ingroup member
        ani = min(v[0] for v in vals)
        af = min(v[1] for v in vals)
        if not (args.min_ani <= ani <= args.max_ani):
            continue
        if af < args.min_af:
            continue
        cands.append({"sample": t, "ani": ani, "af": af})
    # rank by aligned fraction first: placement reliability matters more than
    # raw identity for a module whose whole job is mapping 500 bp anchors
    cands.sort(key=lambda c: (-c["af"], -c["ani"]))
    log("%d candidates inside the band (ANI %.1f-%.1f, AF>=%.0f)"
        % (len(cands), args.min_ani, args.max_ani, args.min_af))
    if not cands:
        log("NO CANDIDATES -- widen --min-ani or lower --min-af")
        open(os.path.join(args.outdir, "outgroups.manifest"), "w").close()
        return 1

    # independence: pairwise ANI among the top candidates
    pool = {sample_name(p): p for p in
            (l.strip() for l in open(args.pool_list) if l.strip())}
    top = cands[: max(args.n_candidates, args.n_outgroups * 8)]
    tl = os.path.join(args.outdir, "top_candidates.txt")
    with open(tl, "w") as fh:
        for c in top:
            if c["sample"] in pool:
                fh.write(pool[c["sample"]] + "\n")
    d2 = os.path.join(args.outdir, "candidate_pairwise.tsv")
    if not (args.reuse and os.path.exists(d2)):
        log("skani: candidate independence ...")
        skani_dist(tl, tl, d2, args.threads, args.skani, 50.0)
    pw = read_dist(d2)

    chosen = []
    for c in top:
        if len(chosen) >= args.n_outgroups:
            break
        dup = False
        for s in chosen:
            v = pw.get((c["sample"], s["sample"])) or pw.get((s["sample"], c["sample"]))
            if v and v[0] >= args.max_outgroup_ani:
                dup = True
                break
        if not dup:
            chosen.append(c)

    with open(os.path.join(args.outdir, "outgroups.tsv"), "w") as fh:
        fh.write("sample\tmin_ani_to_ingroup\tmin_aligned_fraction\trank\tpath\n")
        for i, c in enumerate(chosen):
            fh.write("%s\t%.3f\t%.3f\t%d\t%s\n"
                     % (c["sample"], c["ani"], c["af"], i + 1,
                        pool.get(c["sample"], ".")))
    with open(os.path.join(args.outdir, "outgroups.manifest"), "w") as fh:
        for c in chosen:
            if c["sample"] in pool:
                fh.write(pool[c["sample"]] + "\n")

    log("=" * 70)
    log("SELECTED %d independent outgroups" % len(chosen))
    for c in chosen:
        log("  %-22s ANI=%.2f  AF=%.1f" % (c["sample"], c["ani"], c["af"]))
    if len(chosen) < 2:
        log("  WARNING: fewer than 2 independent outgroups -- polarity at this "
            "cluster will rest on a single assembly and must stay low-confidence")
    log("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
