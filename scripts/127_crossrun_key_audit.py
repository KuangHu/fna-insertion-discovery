#!/usr/bin/env python3
"""Are the cross-run "new" events new SITES, or the same site keyed differently?

An event key is (canonical insert, canonical +/-100bp target context). The
context window is centred on the insertion offset in the SHORT allele, and that
offset is derived from the anchor pair the panel happened to supply. Two panels
with different membership can place the anchors differently, shift the window,
and key one physical insertion twice.

That would inflate `merged - first round` without any new biology, and it is
exactly the failure mode to rule out before quoting a gain.

For every event found ONLY in run B, ask what its insert key does in run A:

  novel_insert        the insert key does not occur in run A at all.
                      A genuinely new element sequence.

  same_insert_new_ctx the insert key IS in run A, but with a different context
                      key. Either the same element at a genuinely different
                      target site, or the SAME site with a shifted window --
                      these two are separated below by comparing the contexts.

This script stops at those two categories -- it only has the hashed keys, and
a hash cannot say whether two windows overlap. Separating "same element, new
site" from "same site, shifted window" needs coordinates, and is done by
128_crossrun_site_audit.py using the junction positions already recorded in
loci.tsv.
"""
import argparse, collections, csv, sys
csv.field_size_limit(sys.maxsize)


def load(p):
    rows = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--run-a-prefix", required=True, help="e.g. full_ecoli")
    ap.add_argument("--run-b-prefix", required=True, help="e.g. deep_ecoli")
    a = ap.parse_args()

    rows = load(a.merged)
    ev = collections.defaultdict(list)
    for r in rows:
        ev[r["event_id"]].append(r)

    # insert keys present in run A, and the context keys they appear with
    a_ctx = collections.defaultdict(set)
    for r in rows:
        if r["locus_id"].split("/", 1)[0] == a.run_a_prefix:
            a_ctx[r["insert_key"]].add(r["context_key"])

    out = collections.Counter()
    for e, rs in ev.items():
        runs = {r["locus_id"].split("/", 1)[0] for r in rs}
        if runs != {a.run_b_prefix}:
            continue
        out["b_only_events"] += 1
        ik = rs[0]["insert_key"]
        if ik not in a_ctx:
            out["novel_insert"] += 1
        elif rs[0]["context_key"] in a_ctx[ik]:
            # should not happen: same insert AND context would have merged
            out["IMPOSSIBLE_same_both"] += 1
        else:
            out["same_insert_new_ctx"] += 1

    n = max(1, out["b_only_events"])
    print("events only in %s : %d" % (a.run_b_prefix, out["b_only_events"]))
    print("  novel insert key            %6d  (%.1f%%)"
          % (out["novel_insert"], 100.0 * out["novel_insert"] / n))
    print("  insert seen in A, new ctx   %6d  (%.1f%%)"
          % (out["same_insert_new_ctx"], 100.0 * out["same_insert_new_ctx"] / n))
    if out["IMPOSSIBLE_same_both"]:
        print("  ** %d events share BOTH keys with run A -- the merge is broken **"
              % out["IMPOSSIBLE_same_both"])
    print()
    print("  A high `new ctx` share means the gain rests on the context window,")
    print("  which is anchor-placement dependent and therefore panel dependent.")
    print("  A high `novel insert` share means the gain is new element sequence,")
    print("  which no window shift can manufacture.")


if __name__ == "__main__":
    sys.exit(main())
