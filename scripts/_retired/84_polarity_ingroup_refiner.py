#!/usr/bin/env python3
"""Layer 3a.5 -- carve a real clade out of a discovery cluster, for polarity only.

Two different jobs need two different groupings, and conflating them is what
blocked Layer 3c:

    Layer 1/2 need genomes similar enough to compare homologous loci.
      Single-linkage chaining is entirely correct for that.
    Layer 3 needs a set with a definable MRCA.
      Chaining is fatal for that.

Measured on cl0000: ingroup pairwise ANI spans 98.92-99.94 while the selected
outgroups reach 98.97 -- the outgroup ceiling sits ABOVE the ingroup floor, so
some outgroups are closer to the ingroup than ingroup members are to each other.
Local trees interleave them and monophyly fails by construction, not by
recombination. 2/30 cl0000 pairs are below the nominal 99.0 clustering
threshold, i.e. stage 10 chained them in.

So this module does NOT touch stage 10. It reads a discovery cluster and emits a
polarity-focal clade defined by COMPLETE linkage (all pairs, no chaining), plus
an explicit separation gap from the nearest external genome:

    ani_gap = min(pairwise ANI inside the clade)
            - max(ANI of any external genome to the clade)

The gap is the real requirement. --min-ani-floor and --max-outgroup-ani are
starting points, not biology: the numbers that matter are reported per clade
(`ingroup_ani_floor`, `closest_outgroup_ani`, `ani_gap`) so the threshold can be
revisited against evidence rather than trusted.

SIZE IS NOT THE OBJECTIVE. A 3-genome clade with a high floor and a clean gap
infers ancestry better than a 6-genome chained set with 5 confusing near
neighbours, so candidates are ranked by gap first and size only as a tiebreak.

Genomes dropped from the focal clade are NOT discarded. They stay valid Layer 2
allele observations and are written to `near_external.manifest`, eligible to act
as extra external lineages in Layer 3c wherever the local genealogy places them
outside the focal clade.

This does not replace the local-tree monophyly gate. Recombination can still
break monophyly at an individual locus even inside a clean clade, so Layer 3c
must keep requiring `local_ingroup_monophyly` and returning unresolved without
it.

    84_polarity_ingroup_refiner.py --cluster cl0000.manifest \\
        --pairwise stage10/skani_dist.tsv --external outgroup/cl0000/ingroup_vs_pool.tsv \\
        --outdir polarity_ingroup/cl0000
"""
import argparse
import csv
import itertools
import os
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


def read_pairs(path):
    """(a, b) -> ANI, symmetric."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, t = sample_name(r["Query_file"]), sample_name(r["Ref_file"])
            try:
                a = float(r["ANI"])
            except (ValueError, KeyError):
                continue
            out[(q, t)] = a
            out.setdefault((t, q), a)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cluster", required=True, help="discovery cluster manifest")
    ap.add_argument("--pairwise", required=True,
                    help="skani dist covering all pairs INSIDE the cluster")
    ap.add_argument("--external", required=True,
                    help="skani dist of cluster members vs the genome pool "
                         "(context only: reports how densely the pool samples "
                         "around the clade)")
    ap.add_argument("--outgroups", default=None,
                    help="outgroups.tsv from 80_outgroup_selector. The gap is "
                         "measured against THESE, because these are the taxa "
                         "that actually enter the local tree. A 99.98%% relative "
                         "sitting unused in a 13k-genome pool cannot break "
                         "monophyly among taxa that are not in the alignment.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-ani-floor", type=float, default=99.5,
                    help="all-pairs ANI floor inside the focal clade")
    ap.add_argument("--min-gap", type=float, default=0.3,
                    help="required separation between the clade floor and the "
                         "closest external genome; THIS is the real criterion")
    ap.add_argument("--min-size", type=int, default=3)
    ap.add_argument("--relax-step", type=float, default=0.25,
                    help="if nothing qualifies, lower the floor by this and retry")
    ap.add_argument("--min-floor-allowed", type=float, default=99.0)
    ap.add_argument("--dup-ani", type=float, default=99.99,
                    help="external genomes at/above this ANI to a cluster "
                         "member are REDEPOSITS, not lineages. NCBI carries "
                         "every assembly twice (GCA_x and GCF_x are the same "
                         "sequence at 100%%), and strains are re-submitted under "
                         "new accessions at ~99.99%%. Counting those as the "
                         "'closest external genome' makes every separation gap "
                         "negative and is measurement error, not biology.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    paths = {sample_name(p): p for p in
             (l.strip() for l in open(args.cluster) if l.strip())}
    members = sorted(paths)
    pw = read_pairs(args.pairwise)
    log("cluster: %d genomes" % len(members))

    # closest external genome to each member (excluding cluster members)
    def core(n):
        """GCA_000258145.1 -> 000258145 ; GCA and GCF share it."""
        p_ = n.split("_")
        return p_[1].split(".")[0] if len(p_) > 1 else n

    member_cores = {core(m) for m in members}
    ext = defaultdict(float)
    n_dup = 0
    ext_best = {}
    with open(args.external) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, t = sample_name(r["Query_file"]), sample_name(r["Ref_file"])
            if q not in paths or t in paths:
                continue
            try:
                a = float(r["ANI"])
            except (ValueError, KeyError):
                continue
            if core(t) in member_cores or a >= args.dup_ani:
                n_dup += 1
                continue                      # redeposit, not a lineage
            if a > ext[q]:
                ext[q] = a
                ext_best[q] = t
    log("excluded %d redeposit/duplicate comparisons (>=%.2f ANI or shared "
        "accession core)" % (n_dup, args.dup_ani))

    def floor_of(sub):
        vs = [pw.get((a, b)) for a, b in itertools.combinations(sub, 2)]
        vs = [v for v in vs if v is not None]
        return min(vs) if vs else None

    # The gap that matters is against the taxa that will actually be in the
    # tree. Measured against the whole pool it is always negative -- the pool
    # is dense enough that some E. coli sits at 99.98% of any cluster -- but
    # that genome is never placed in the alignment, so it cannot interleave.
    sel = {}
    if args.outgroups and os.path.exists(args.outgroups):
        with open(args.outgroups) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    sel[r["sample"]] = float(r["min_ani_to_ingroup"])
                except (ValueError, KeyError):
                    continue
    sel_max = max(sel.values()) if sel else None
    pool_closest = max(ext.values()) if ext else 0.0
    if sel_max is not None:
        log("selected outgroups: %d, closest at %.3f ANI (pool closest %.3f, "
            "not used)" % (len(sel), sel_max, pool_closest))

    best, floor = None, args.min_ani_floor
    while floor >= args.min_floor_allowed:
        cands = []
        for k in range(len(members), args.min_size - 1, -1):
            for sub in itertools.combinations(members, k):
                f = floor_of(sub)
                if f is None or f < floor:
                    continue                      # COMPLETE linkage: all pairs
                closest = (sel_max if sel_max is not None
                           else max((ext.get(m, 0.0) for m in sub), default=0.0))
                gap = f - closest
                cands.append({"members": list(sub), "floor": f,
                              "closest": closest, "gap": gap})
        # gap first, size only as a tiebreak: a clean 3-genome clade beats a
        # chained 6-genome set
        cands.sort(key=lambda c: (-c["gap"], -len(c["members"])))
        ok = [c for c in cands if c["gap"] >= args.min_gap]
        if ok:
            best = ok[0]
            break
        if cands and best is None:
            best = cands[0]                       # remember the best seen
        floor -= args.relax_step

    with open(os.path.join(args.outdir, "refinement.tsv"), "w") as fh:
        fh.write("discovery_cluster_size\tfocal_clade_size\tingroup_ani_floor\t"
                 "closest_outgroup_ani\tani_gap\tgap_requirement_met\t"
                 "pool_closest_ani\tfocal_members\texcluded_members\n")
        if best:
            excl = [m for m in members if m not in best["members"]]
            fh.write("%d\t%d\t%.3f\t%.3f\t%.3f\t%s\t%.3f\t%s\t%s\n"
                     % (len(members), len(best["members"]), best["floor"],
                        best["closest"], best["gap"],
                        "yes" if best["gap"] >= args.min_gap else "NO",
                        pool_closest,
                        ",".join(best["members"]), ",".join(excl) or "."))
        else:
            fh.write("%d\t0\tNA\tNA\tNA\tNO\tNA\t.\t%s\n"
                     % (len(members), ",".join(members)))

    if not best:
        log("NO focal clade found at any floor >= %.2f" % args.min_floor_allowed)
        for n in ("focal_clade.manifest", "near_external.manifest"):
            open(os.path.join(args.outdir, n), "w").close()
        return 1

    excl = [m for m in members if m not in best["members"]]
    with open(os.path.join(args.outdir, "focal_clade.manifest"), "w") as fh:
        for m in best["members"]:
            fh.write(paths[m] + "\n")
    # dropped genomes stay valid Layer 2 observations and may serve as extra
    # external lineages wherever the local tree puts them outside the clade
    with open(os.path.join(args.outdir, "near_external.manifest"), "w") as fh:
        for m in excl:
            fh.write(paths[m] + "\n")

    log("=" * 72)
    log("POLARITY INGROUP REFINEMENT")
    log("  discovery cluster      : %d genomes" % len(members))
    log("  focal clade            : %d genomes  %s"
        % (len(best["members"]), ",".join(best["members"])))
    log("  ingroup_ani_floor      : %.3f  (all pairs, complete linkage)"
        % best["floor"])
    log("  closest_outgroup_ani   : %.3f  (selected outgroups)" % best["closest"])
    log("  pool_closest_ani       : %.3f  (context only, not in the tree)"
        % pool_closest)
    log("  ani_gap                : %+.3f  %s"
        % (best["gap"], "OK" if best["gap"] >= args.min_gap
           else "BELOW REQUIREMENT (%.2f)" % args.min_gap))
    if excl:
        log("  moved to near_external : %s" % ",".join(excl))
    if best["gap"] < args.min_gap:
        log("")
        log("  WARNING: no positive separation was achievable. Polarity from "
            "this cluster will remain unreliable however good the trees look.")
    log("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
