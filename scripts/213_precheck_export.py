#!/usr/bin/env python3
"""Preconditions the v1 record format needs from v2, measured not assumed.

Three things must hold before an export can claim "every sequence is replayable
from its declared coordinates":

  P1  The published reconstructed sequence is EXACTLY
      carrier[s-up_len : s] + carrier[e : e+dn_len] at the coordinates the row
      records. If it is not, the row's declared source interval is wrong and no
      amount of re-export fixes it.

  P2  chosen_shift is REPORTED, not APPLIED. v2 stores the original interval in
      insert_start/insert_end and publishes the sequence cut there; a nonzero
      chosen_shift names a verified-equivalent alternative, it does not move the
      exported bases. The export must therefore set chosen_joint_shift_bp = 0
      and carry the nonzero value as range information, or the per-flank slice
      check in section 8C would fail even though the joined product is right.

  P3  Flank lengths are clipped at contig ends, so up_len/dn_len must be read
      per row and never assumed to be 2000.
"""
import argparse
import collections
import random
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"
V2 = W + "/inferred_empty_v2"


def fa(p):
    o, k, b = {}, None, []
    with open(p) as fh:
        for l in fh:
            if l[0] == ">":
                if k is not None:
                    o[k] = "".join(b).upper()
                k = l[1:].split()[0]
                b = []
            else:
                b.append(l.strip())
    if k is not None:
        o[k] = "".join(b).upper()
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=20261108)
    a = ap.parse_args()

    man = {}
    with open(SRC + "/genome_map.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                man[f[0]] = f[1]

    rows = []
    with open(V2 + "/occurrences.tsv") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["official_fasta"]] != "yes":
                continue
            rows.append(tuple(f[ix[c]] for c in (
                "occurrence_id", "carrier_accession", "carrier_contig",
                "insert_start", "insert_end", "carrier_up_len",
                "carrier_dn_len", "chosen_shift")))
    rnd = random.Random(a.seed)
    shifted = [r for r in rows if int(r[7]) != 0]
    sample = rnd.sample(rows, min(a.n, len(rows))) + shifted
    want = {r[0] for r in sample}

    pub, cur, buf = {}, None, []
    with open(V2 + "/inferred_empty.fna") as fh:
        for l in fh:
            if l[0] == ">":
                if cur in want:
                    pub[cur] = "".join(buf)
                cur = l[1:].split()[0]
                buf = []
            elif cur in want:
                buf.append(l.strip())
    if cur in want:
        pub[cur] = "".join(buf)

    by = collections.defaultdict(list)
    for r in sample:
        by[r[1]].append(r)

    res = collections.Counter()
    bad = []
    shift_applied = collections.Counter()
    for acc, rs in by.items():
        p = man.get(acc)
        if not p:
            res["carrier_not_in_manifest"] += len(rs)
            continue
        g = fa(p)
        for occ, _, ctg, s, e, ul, dl, sh in rs:
            S = g.get(ctg)
            if S is None:
                res["contig_missing"] += 1
                continue
            s, e, ul, dl, sh = int(s), int(e), int(ul), int(dl), int(sh)
            P = pub.get(occ)
            if P is None:
                res["no_published_record"] += 1
                continue
            at0 = S[max(0, s - ul):s] + S[e:e + dl]
            if at0 == P:
                res["P1_slice_at_recorded_interval_matches"] += 1
                if sh != 0:
                    shift_applied["nonzero_shift_but_sequence_cut_at_shift_0"] += 1
            else:
                k = sh
                atk = S[max(0, s + k - ul):s + k] + S[e + k:e + k + dl]
                if k and atk == P:
                    res["P1_matches_only_with_shift_applied"] += 1
                    shift_applied["nonzero_shift_IS_applied"] += 1
                else:
                    res["P1_MISMATCH"] += 1
                    bad.append((occ, len(P), len(at0)))
            if len(S) >= s and (ul != min(2000, s) or
                                dl != min(2000, max(0, len(S) - e))):
                res["P3_flank_len_not_the_clipped_default"] += 1

    print("sampled official rows: %d  (plus every nonzero-shift row: %d)"
          % (len(sample) - len(shifted), len(shifted)))
    for k, v in sorted(res.items()):
        print("   %-46s %7d" % (k, v))
    print("   -- shift semantics --")
    for k, v in sorted(shift_applied.items()):
        print("   %-46s %7d" % (k, v))
    for b in bad[:10]:
        print("   MISMATCH %s pub_len=%d rebuilt_len=%d" % b)
    if res["P1_MISMATCH"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
