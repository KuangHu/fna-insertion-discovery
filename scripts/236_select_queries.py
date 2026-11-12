#!/usr/bin/env python3
"""Choose the DNA query set: each bag's original representative, plus the
variants that representative genuinely fails to cover.

THE RULE (unchanged from the design the user corrected into place):
  * Every bag KEEPS ITS ORIGINAL REPRESENTATIVE -- `cluster_rep` from
    insertion_cds_strict.tsv, the same sequence stage 232 searched with.
    Nothing replaces it, not even a longer member.
  * A member is absorbed only if an already-kept representative of that bag
    covers it TWO-SIDED: >=MINID identity AND >=MINCOV of the MEMBER AND
    >=MINCOV of the REPRESENTATIVE, measured on merged HSP blocks on each
    sequence separately.  One-sided (shorter-sequence) coverage is what made
    the first attempt wrong: a 1.3 kb IS inside a 5 kb Tn3 covers 100% of the
    IS and ~26% of the Tn3, so single linkage on the shorter side fused two
    distinct insertion structures and the IS stopped being searched at all.
    Under the two-sided test that member is KEPT.
  * Anything not absorbed is promoted to an additional representative, so a
    bag whose DNA is heterogeneous is searched with several queries.

WHY THIS IS STREAMED AND NOT HELD IN A DICT.  Bags are near-complete
similarity graphs: measured pair density is 0.89-1.00 of n^2, so the 6,144
bags hold ~1e8 aligned pairs and the largest single bag ~4.9e7.  A dict keyed
by pair is hundreds of GB and dies.  Blast groups its output by query, so this
indexes each bag's alignment file by query id once and then reads one member's
rows at a time: peak memory is one member's block, and BOTH coverages come out
of that single block because the format carries qstart/qend and sstart/send
together.  Reading the member's block alone is not enough, for two separate reasons
measured in this dataset.  In 1 pair blast reported the alignment only as
(representative -> member), so the member's block held nothing.  In 2 further
pairs it reported BOTH directions but with DIFFERENT HSP sets, the member's
own row being the weaker of the two -- so preferring the member's row whenever
it exists still discards the better evidence.  A direction that is missing or
weaker is a detection artifact, not evidence of non-containment.  So each
representative's own block is cached when it is promoted, and the pair is
judged on whichever direction demonstrates containment best.  A representative
can only absorb members considered after it was promoted, which keeps the
result independent of the order blast happened to report things in.

COUNTS THIS FILE DOES NOT MAKE.  A representative is a SEARCH QUERY.  It is
not a copy, not an insertion site, and not a transposition event.  Absorbing
one member into another says their DNA is mutually contained at MINCOV; it
says nothing about how many times the element moved.
"""
import argparse, collections, csv, json, os


def merged_len(iv):
    iv.sort(); tot = 0; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    return tot + ce - cs + 1


def index_by_query(path):
    """qseqid -> list of (byte_start, byte_end). Blast groups rows by query, so
    the common case is one contiguous range; ranges are merged when adjacent so
    a file that is not grouped still indexes correctly rather than silently
    losing rows."""
    idx = collections.defaultdict(list); off = 0
    with open(path, "rb") as f:
        for line in f:
            n = len(line); q = line.split(b"\t", 1)[0]
            r = idx[q]
            if r and r[-1][1] == off: r[-1][1] = off + n
            else: r.append([off, off + n])
            off += n
    return idx


def rows_for(fh, ranges):
    for s, e in ranges:
        fh.seek(s)
        for line in fh.read(e - s).decode().splitlines():
            p = line.split("\t")
            if len(p) >= 10: yield p


def coverage_against(fh, ranges, minid):
    """-> {subject: (cov_of_query, cov_of_subject)} from one query's own block."""
    qb = collections.defaultdict(list); sb = collections.defaultdict(list)
    ql = {}; sl = {}
    for p in rows_for(fh, ranges):
        q, s, pid = p[0], p[1], float(p[2])
        if q == s or pid < minid: continue
        qs, qe, ss, se = int(p[4]), int(p[5]), int(p[6]), int(p[7])
        if qe < qs: qs, qe = qe, qs
        if se < ss: ss, se = se, ss
        qb[s].append((qs, qe)); sb[s].append((ss, se))
        ql[s] = int(p[8]); sl[s] = int(p[9])
    return {s: (merged_len(qb[s]) / ql[s], merged_len(sb[s]) / sl[s]) for s in qb}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bagaln", required=True, help="dir of <bag>.tsv alignments")
    ap.add_argument("--bagfa", required=True, help="dir of <bag>.fna")
    ap.add_argument("--cds", required=True, help="insertion_cds_strict.tsv")
    ap.add_argument("--occbags", required=True, help="strict_occurrence_bags.tsv")
    ap.add_argument("--fna", required=True, help="dedup_insertions.fna")
    ap.add_argument("--out", required=True, help="output dir")
    ap.add_argument("--minid", type=float, default=90.0)
    ap.add_argument("--mincov", type=float, default=0.90)
    a = ap.parse_args()

    bagkeys = collections.defaultdict(set)
    with open(a.occbags) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bagkeys[r["cds_cluster_id"]].add(r["seq_key"])

    orig = {}
    with open(a.cds) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b = r.get("cds_cluster_id")
            if b and b not in orig:
                orig[b] = r.get("cluster_rep") or r.get("seq_key")

    lens = {}
    k = None; n = 0
    with open(a.fna) as f:
        for line in f:
            if line[0] == ">":
                if k: lens[k] = n
                k = line[1:].split()[0].replace("SEQ|", ""); n = 0
            else: n += len(line.strip())
        if k: lens[k] = n

    repaired = []
    for b in list(bagkeys):
        bagkeys[b] = {x for x in bagkeys[b] if x in lens}
        if not bagkeys[b]: del bagkeys[b]; continue
        if b not in orig or orig[b] not in bagkeys[b]:
            repaired.append((b, orig.get(b, "<absent>"), sorted(bagkeys[b])[0]))
            orig[b] = sorted(bagkeys[b])[0]

    reps = {}; decisions = []
    for b in sorted(bagkeys):
        keys = bagkeys[b]; o = orig[b]; reps[b] = [o]
        p = os.path.join(a.bagaln, b + ".tsv")
        if len(keys) == 1 or not os.path.exists(p) or os.path.getsize(p) == 0:
            for m in sorted(keys - {o}):
                reps[b].append(m)
                decisions.append((b, m, "kept_as_variant", "", 0.0, 0.0))
            continue
        idx = index_by_query(p)
        with open(p, "rb") as fh:
            # a representative's own block, cached when it is promoted, so the
            # pair is seen even if blast reported only that direction
            fwd = {o: coverage_against(fh, idx.get(o.encode(), []), a.minid)}
            # same deterministic order as the rule: longest first, then key
            for m in sorted(keys - {o}, key=lambda x: (-lens[x], x)):
                cov = coverage_against(fh, idx.get(m.encode(), []), a.minid)
                hit = None
                for r in reps[b]:
                    # blast may report BOTH directions with different HSP sets,
                    # so take whichever demonstrates containment best rather
                    # than preferring the member's own row when it exists
                    cands = []
                    if cov.get(r) is not None: cands.append(cov[r])
                    rc = fwd[r].get(m)
                    if rc is not None: cands.append((rc[1], rc[0]))  # (cov_member, cov_rep)
                    if not cands: continue
                    c = max(cands, key=lambda x: min(x[0], x[1]))
                    if c[0] >= a.mincov and c[1] >= a.mincov: hit = (r, c); break
                if hit:
                    decisions.append((b, m, "absorbed", hit[0],
                                      round(hit[1][0], 4), round(hit[1][1], 4)))
                else:
                    reps[b].append(m)
                    fwd[m] = coverage_against(fh, idx.get(m.encode(), []), a.minid)
                    decisions.append((b, m, "kept_as_variant", "", 0.0, 0.0))

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "query_decisions.tsv"), "w") as f:
        f.write("bag_id\tseq_key\tdecision\tabsorbed_by\tcov_of_member\tcov_of_representative\n")
        for d in decisions: f.write("\t".join(str(x) for x in d) + "\n")
    with open(os.path.join(a.out, "query_repaired_reps.tsv"), "w") as f:
        f.write("bag_id\toriginal_rep_as_recorded\tsubstitute\n")
        for d in repaired: f.write("\t".join(d) + "\n")

    seqs = {}
    want = {x for v in reps.values() for x in v}
    k = None; buf = []
    with open(a.fna) as f:
        for line in f:
            if line[0] == ">":
                if k and k in want: seqs[k] = "".join(buf)
                k = line[1:].split()[0].replace("SEQ|", ""); buf = []
            else: buf.append(line.strip())
        if k and k in want: seqs[k] = "".join(buf)

    rows = [("Q_" + x, b, x, i == 0) for b in sorted(reps) for i, x in enumerate(reps[b])]
    seen = set()
    with open(os.path.join(a.out, "queries.fna"), "w") as f:
        for qid, b, x, _ in rows:
            if x in seen: continue
            seen.add(x); f.write(">%s\n%s\n" % (qid, seqs[x]))
    with open(os.path.join(a.out, "queries.tsv"), "w") as f:
        f.write("query_id\tbag_id\tseq_key\tis_bag_original_rep\tlen_nt\n")
        for qid, b, x, io in rows:
            f.write("%s\t%s\t%s\t%d\t%d\n" % (qid, b, x, int(io), lens[x]))

    nv = {b: len(v) for b, v in reps.items()}
    summ = dict(
        bags=len(reps), query_records_bag_x_rep=len(rows),
        distinct_query_sequences=len(seen),
        total_query_bp=sum(lens[x] for x in seen),
        bags_with_one_representative=sum(1 for c in nv.values() if c == 1),
        bags_with_multiple_representatives=sum(1 for c in nv.values() if c > 1),
        max_representatives_for_one_bag=max(nv.values()),
        members_absorbed_two_sided=sum(1 for d in decisions if d[2] == "absorbed"),
        members_kept_as_distinct_variants=sum(1 for d in decisions if d[2] == "kept_as_variant"),
        original_reps_repaired=len(repaired))
    with open(os.path.join(a.out, "queries.json"), "w") as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
