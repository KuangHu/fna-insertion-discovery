#!/usr/bin/env python3
"""Same element at a NEW SITE, or the SAME SITE keyed twice?

127_ split the cross-run gain into events whose insert sequence is novel (which
no window shift can manufacture) and events whose insert was already seen in
run A under a different context key. This resolves that second group.

The context key is a hash, so it cannot be compared for overlap. But every
decomposed locus records where its junction sits in a real genome:

    pre_event_genome, pre_event_contig, junction_L_abs, junction_R_abs

If a run-B locus and a run-A locus carrying the SAME insert key land on the
same contig of the same genome within TOL bp, they are the same physical site
described twice, and the run-B event is an artefact of the key rather than a
discovery. If they land far apart, or in different genomes with no shared
frame, the element is genuinely at another site -- which is what a mobile
element does, and is a real finding.

  same_site        a run-A locus with the same insert key is within TOL bp
                   on the same contig of the same genome  -> INFLATION
  different_site   a shared frame exists and the junctions are far apart
                   -> real, the same element elsewhere
  no_shared_frame  no run-A locus with that insert key resolved a junction in
                   a genome this run-B locus also resolved -> undecidable here

Only events already classified `insert seen in A, new ctx` are examined;
novel-insert events cannot be the same site by construction.
"""
import argparse, collections, csv, glob, os, sys
csv.field_size_limit(sys.maxsize)


def junctions(run_dir, run_prefix):
    """locus_id -> (genome, contig, jL, jR) for every placed locus."""
    out = {}
    for d in sorted(glob.glob(os.path.join(run_dir, "l2", "*/"))):
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        pl = "%s/%s" % (run_prefix, os.path.basename(d.rstrip("/")))
        for r in csv.DictReader(open(p), delimiter="\t"):
            g, c = r.get("pre_event_genome", "."), r.get("pre_event_contig", ".")
            if g in (".", "") or c in (".", ""):
                continue
            try:
                jl, jr = int(r["junction_L_abs"]), int(r["junction_R_abs"])
            except (ValueError, KeyError):
                continue
            out["%s|%s" % (pl, r["locus_id"])] = (g, c, jl, jr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--run-a", required=True)
    ap.add_argument("--run-b", required=True)
    ap.add_argument("--tol", type=int, default=500)
    a = ap.parse_args()
    pa, pb = os.path.basename(a.run_a.rstrip("/")), os.path.basename(a.run_b.rstrip("/"))

    rows = list(csv.DictReader(open(a.merged), delimiter="\t"))
    ev = collections.defaultdict(list)
    for r in rows:
        ev[r["event_id"]].append(r)

    # run-A loci indexed by insert key
    a_by_ik = collections.defaultdict(list)
    a_ctx = collections.defaultdict(set)
    for r in rows:
        if r["locus_id"].split("/", 1)[0] == pa:
            a_by_ik[r["insert_key"]].append(r["locus_id"])
            a_ctx[r["insert_key"]].add(r["context_key"])

    ja, jb = junctions(a.run_a, pa), junctions(a.run_b, pb)
    out = collections.Counter()
    for e, rs in ev.items():
        if {r["locus_id"].split("/", 1)[0] for r in rs} != {pb}:
            continue
        ik = rs[0]["insert_key"]
        if ik not in a_ctx or rs[0]["context_key"] in a_ctx[ik]:
            continue
        out["ambiguous_events"] += 1
        verdict = "no_shared_frame"
        for r in rs:
            mine = jb.get(r["locus_id"])
            if not mine:
                continue
            for al in a_by_ik[ik]:
                theirs = ja.get(al)
                if not theirs:
                    continue
                if mine[0] == theirs[0] and mine[1] == theirs[1]:
                    if abs(mine[2] - theirs[2]) <= a.tol:
                        verdict = "same_site"
                        break
                    verdict = "different_site"
            if verdict == "same_site":
                break
        out[verdict] += 1

    n = max(1, out["ambiguous_events"])
    print("ambiguous events (insert in A, new context) : %d" % out["ambiguous_events"])
    for k in ("same_site", "different_site", "no_shared_frame"):
        print("  %-18s %6d  (%.1f%%)" % (k, out[k], 100.0 * out[k] / n))
    print()
    print("  same_site is the INFLATION: one physical insertion keyed twice")
    print("  because the anchor pair moved the context window.")


if __name__ == "__main__":
    sys.exit(main())
