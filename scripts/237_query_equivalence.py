#!/usr/bin/env python3
"""Record which DNA queries from DIFFERENT bags are the same element.

This is deliberately a RECORD, not an action.  Bags are cut on protein
similarity and the search runs on DNA, so two bags can carry the same element:
in the 181-reference pilot, every one of BAGS000173's 2,861 loci was also hit
by BAGS002634.  The tempting fix -- delete one of the queries -- is wrong,
because a query dropped up front can never be recovered, and two bags that
share one element may still differ elsewhere.  So both queries keep searching
and the equivalence is applied later at the LOCUS level, where collapsing two
names for one place costs no sensitivity.

Equivalence uses the same two-sided test as the representative choice:
>=MINID identity and >=MINCOV of BOTH sequences on merged HSP blocks.  A query
contained inside a larger one is NOT equivalent to it and is not recorded as
such.

Streamed per query for the same reason as 236: the query set is highly
self-similar, so a dict over all pairs is not affordable.
"""
import argparse, collections, csv, importlib.util, json, os

_s = importlib.util.spec_from_file_location(
    "sel", os.path.join(os.path.dirname(os.path.abspath(__file__)), "236_select_queries.py"))
sel = importlib.util.module_from_spec(_s); _s.loader.exec_module(sel)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True, help="queries_selfaln.tsv")
    ap.add_argument("--queries", required=True, help="queries.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--minid", type=float, default=90.0)
    ap.add_argument("--mincov", type=float, default=0.90)
    a = ap.parse_args()

    bag_of = {}
    with open(a.queries) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bag_of.setdefault(r["query_id"], r["bag_id"])

    idx = sel.index_by_query(a.aln)
    pairs = {}
    with open(a.aln, "rb") as fh:
        for q in idx:
            qs = q.decode()
            for s, (ca, cb) in sel.coverage_against(fh, idx[q], a.minid).items():
                if ca >= a.mincov and cb >= a.mincov:
                    k = (qs, s) if qs < s else (s, qs)
                    v = (ca, cb) if qs < s else (cb, ca)
                    old = pairs.get(k)
                    if old is None or (v[0] + v[1]) > (old[0] + old[1]): pairs[k] = v

    cross = 0
    with open(a.out, "w") as f:
        f.write("query_a\tquery_b\tbag_a\tbag_b\tsame_bag\tcov_a\tcov_b\n")
        for (x, y), (ca, cb) in sorted(pairs.items()):
            ba, bb = bag_of.get(x, ""), bag_of.get(y, "")
            same = int(ba == bb and ba != "")
            cross += 0 if same else 1
            f.write("%s\t%s\t%s\t%s\t%d\t%.4f\t%.4f\n" % (x, y, ba, bb, same, ca, cb))

    comp = {}
    def find(z):
        while comp.get(z, z) != z: comp[z] = comp.get(comp[z], comp[z]); z = comp[z]
        return z
    for (x, y) in pairs:
        rx, ry = find(x), find(y)
        if rx != ry: comp[max(rx, ry)] = min(rx, ry)
    groups = collections.defaultdict(set)
    for q in bag_of: groups[find(q)].add(bag_of[q])
    multibag = sum(1 for g in groups.values() if len(g) > 1)

    summ = dict(queries=len(bag_of), equivalent_pairs=len(pairs),
                equivalent_pairs_across_different_bags=cross,
                equivalence_components=len(groups),
                components_spanning_more_than_one_bag=multibag,
                note="recorded only; no query is removed because of this file")
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
