#!/usr/bin/env python3
"""Where does the support for a cross-run NEW event actually come from?

The deep run added 299 events the first round did not have, while introducing
only 8 genomes. That is CONSISTENT with "new comparisons among old genomes"
but does not show it: 8 genomes can carry many alleles nothing else had.

This asks the question the existing output can actually answer. For each new
event, look at which genomes carry its two decomposed alleles and ask whether
the first round already held them:

  common_supportable   BOTH the long (filled) and the short (empty) allele have
                       at least one carrier among the genomes SHARED by the two
                       runs. The difference was therefore visible in genomes the
                       first round already used, and was missed for some other
                       reason -- panel composition, graph construction, or the
                       anchor pair -- which this does NOT separate.

  needs_new_genomes    at least one side is carried ONLY by genomes the deep run
                       introduced. Without those genomes the comparison could
                       not have been made at all.

  unattributable       carrier lists incomplete, or placement ambiguous, so the
                       event cannot be assigned either way.

An event is scored on its BEST locus: if any one of its loci is supportable
from shared genomes, the event is. That rule is GENEROUS to
`common_supportable`, so read the two numbers as bounds in opposite directions:

    common_supportable   an UPPER bound on what was reachable from shared genomes
    needs_new_genomes    a LOWER bound on what required the new genomes

(An earlier version of this docstring argued the opposite and was wrong: a
generous rule inflates its favoured category, it does not protect it.)

What this is NOT: causal attribution. It bounds how much of the gain COULD have
come from genomes already in hand. Showing that a given event was missed
BECAUSE of panel composition rather than graph construction needs a
counterfactual run, not a carrier list.
"""
import argparse, collections, csv, glob, os, sys
csv.field_size_limit(sys.maxsize)


def manifest_genomes(d):
    g = set()
    for p in glob.glob(os.path.join(d, "panels_manifests", "*.manifest")):
        for line in open(p):
            b = os.path.basename(line.strip())
            if b:
                g.add(b[:-4] if b.endswith(".fna") else b)
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-a", required=True, help="first-round run dir")
    ap.add_argument("--run-b", required=True, help="deep run dir")
    ap.add_argument("--merged", required=True, help="merged loci_to_event.tsv")
    a = ap.parse_args()

    ga, gb = manifest_genomes(a.run_a), manifest_genomes(a.run_b)
    common, newg = ga & gb, gb - ga
    print("genomes: run A %d, run B %d, SHARED %d, new in B %d"
          % (len(ga), len(gb), len(common), len(newg)))

    ev_loci = collections.defaultdict(list)
    for r in csv.DictReader(open(a.merged), delimiter="\t"):
        ev_loci[r["event_id"]].append(r["locus_id"])
    ra, rb = os.path.basename(a.run_a.rstrip("/")), os.path.basename(a.run_b.rstrip("/"))
    new_ev = [e for e, ls in ev_loci.items()
              if all(l.split("/", 1)[0] == rb for l in ls)]
    only_a = [e for e, ls in ev_loci.items()
              if all(l.split("/", 1)[0] == ra for l in ls)]
    print("events only in B (the gain)  %d ; only in A  %d ; total %d"
          % (len(new_ev), len(only_a), len(ev_loci)))

    # carriers + placement for every locus of run B
    car, info = collections.defaultdict(dict), {}
    for d in sorted(glob.glob(os.path.join(a.run_b, "l2", "*/"))):
        panel = "%s/%s" % (rb, os.path.basename(d.rstrip("/")))
        lp, apth = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.tsv")
        if not (os.path.exists(lp) and os.path.exists(apth)):
            continue
        for r in csv.DictReader(open(lp), delimiter="\t"):
            info["%s|%s" % (panel, r["locus_id"])] = (
                r["shortest_allele"], r["longest_allele"],
                r.get("placement_status", ""))
        for r in csv.DictReader(open(apth), delimiter="\t"):
            gs = {x for x in (r.get("genomes") or "").split(",") if x}
            car["%s|%s" % (panel, r["locus_id"])][r["allele_id"]] = gs

    def classify(lid):
        meta = info.get(lid)
        if not meta:
            return "unattributable"
        sa, la, st = meta
        cs, cl = car.get(lid, {}).get(sa), car.get(lid, {}).get(la)
        if not cs or not cl:
            return "unattributable"
        if st == "ambiguous":
            return "unattributable"
        if (cs & common) and (cl & common):
            return "common_supportable"
        if (cs and not (cs - newg)) or (cl and not (cl - newg)):
            return "needs_new_genomes"
        return "unattributable"

    RANK = {"common_supportable": 0, "needs_new_genomes": 1, "unattributable": 2}
    out = collections.Counter()
    for e in new_ev:
        out[min((classify(l) for l in ev_loci[e]), key=lambda x: RANK[x])] += 1

    print()
    print("SUPPORT SOURCE of the %d events only run B found" % len(new_ev))
    for k in ("common_supportable", "needs_new_genomes", "unattributable"):
        print("  %-22s %5d  (%.1f%%)"
              % (k, out[k], 100.0 * out[k] / max(1, len(new_ev))))
    print()
    print("  `common_supportable` bounds ABOVE how much of the gain was already")
    print("  reachable from genomes the first round used. It does not say WHY")
    print("  the first round missed it.")


if __name__ == "__main__":
    sys.exit(main())
