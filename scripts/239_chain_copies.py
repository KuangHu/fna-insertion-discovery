#!/usr/bin/env python3
"""Chain blastn HSPs into copies, for a sweep run with MANY QUERIES PER BAG.

233 did this when every bag had exactly one query and the FASTA id WAS the bag
id.  The all-sequence sweep searches with all 63,095 distinct insertion
sequences, so a bag has up to 82 queries and the two identifiers have to be
kept apart: `query_id` is what blast reports, `bag_id` is what the result
belongs to.  Keying by bag id here, as 233 does, would collapse 63,095 queries
onto 6,144 keys and hand every copy the wrong query length.

The chaining and classification are deliberately IDENTICAL to 233, so results
from the two sweeps stay comparable:

  * HSPs on one subject are chained when collinear and separated by less than
    one query length, because a real copy of a 1.3 kb IS can come back as
    three HSPs split by SNP clusters -- while two copies 4 kb apart must not
    become one.
  * full_copy >= 0.90 of the query, partial_copy >= 0.50, else fragment.
  * A copy is called `catalogued` when it reciprocally overlaps a catalogue
    occurrence on the same (assembly, contig) by >= 0.50 both ways.

WHAT A COPY IS NOT.  A copy is not an insertion site and not a transposition
event.  An element present in 10,000 genomes at one ancestral locus is one
locus and zero observed events, not 10,000 of either.  Collapsing copies to
loci happens in the next stage; polarity is not available at all.  Every
column here says `copy`.

CROSS-QUERY REDUNDANCY IS EXPECTED AND IS NOT RESOLVED HERE.  Several queries
of one bag, and queries of different bags that carry the same element, will
all hit the same place.  That is why the sweep is sensitive.  Deduplication is
a locus-level operation and doing it here, per shard, would be wrong: a shard
sees only part of the corpus.
"""
import argparse, collections, csv, json, os, sys


def merged_len(iv):
    iv = sorted(iv); tot = 0; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce: ce = max(ce, e)
        else: tot += ce - cs; cs, ce = s, e
    return tot + ce - cs


def chain(hsps, qlen):
    """hsps: (qs,qe,ss,se,pid,alen,strand) -> [(strand, [hsp,...]), ...]"""
    out = []
    for strand in ("+", "-"):
        h = sorted([x for x in hsps if x[6] == strand], key=lambda x: x[2])
        cur = []
        for x in h:
            if not cur: cur = [x]; continue
            gap = x[2] - max(y[3] for y in cur)
            if gap <= qlen: cur.append(x)
            else: out.append((strand, cur)); cur = [x]
        if cur: out.append((strand, cur))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="blast outfmt 6 with 14 std fields")
    ap.add_argument("--queries", required=True, help="query_id bag_id seq_key len_nt")
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--occurrences", default="", help="catalogue browse tsv (optional)")
    ap.add_argument("--occbags", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--full-qcov", type=float, default=0.90)
    ap.add_argument("--part-qcov", type=float, default=0.50)
    ap.add_argument("--recip", type=float, default=0.50)
    a = ap.parse_args()

    qbag, qlen = {}, {}
    with open(a.queries) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            qbag[r["query_id"]] = r["bag_id"]; qlen[r["query_id"]] = int(r["len_nt"])
    sys.stderr.write("queries %d over %d bags\n" % (len(qbag), len(set(qbag.values()))))

    species = {}
    with open(a.corpus) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            species[r["assembly"]] = r["species"]

    occbag = {}
    if a.occbags and os.path.exists(a.occbags):
        with open(a.occbags) as f:
            rd = csv.DictReader(f, delimiter="\t")
            bcol = ("bag_id" if "bag_id" in rd.fieldnames else
                    "cds_cluster_id" if "cds_cluster_id" in rd.fieldnames else None)
            if not bcol: sys.exit("FATAL: occbags lacks bag_id/cds_cluster_id")
            for r in rd: occbag[r["occurrence_id"]] = r[bcol]

    cat = collections.defaultdict(list)
    if a.occurrences and os.path.exists(a.occurrences):
        with open(a.occurrences) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cat[(r["carrier_assembly_accession_version"],
                     r["carrier_fasta_record_id"])].append(
                    (int(r["insert_start0"]), int(r["insert_end0"]),
                     r["occurrence_id"], occbag.get(r["occurrence_id"], "")))
        for k in cat: cat[k].sort()
    sys.stderr.write("catalogue (assembly,contig) keys %d\n" % len(cat))

    by = collections.defaultdict(list); nh = 0; unknown = set()
    with open(a.hits) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if len(x) < 12: continue
            q, s = x[0], x[1]
            if q not in qlen: unknown.add(q); continue
            qs, qe = int(x[6]), int(x[7]); ss, se = int(x[8]), int(x[9])
            strand = "+" if se >= ss else "-"
            if se < ss: ss, se = se, ss
            by[(q, s)].append((qs, qe, ss, se, float(x[2]), int(x[3]), strand))
            nh += 1
    if unknown:
        sys.exit("FATAL: %d query ids in hits are absent from %s (e.g. %s)"
                 % (len(unknown), a.queries, sorted(unknown)[:3]))
    sys.stderr.write("HSPs %d over %d (query,subject) pairs\n" % (nh, len(by)))

    cols = ["query_id", "bag_id", "assembly", "species", "contig", "s0", "s1", "strand",
            "n_hsp", "qcov", "pident", "copy_kind", "catalogue_status",
            "catalogue_occurrence", "catalogue_bag", "same_bag"]
    n = collections.Counter()
    with open(a.out + ".tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for (q, s), hsps in by.items():
            L = qlen[q]; bag = qbag[q]
            acc, _, contig = s.partition("|")
            for strand, grp in chain(hsps, L):
                qcov = merged_len([(x[0], x[1]) for x in grp]) / L
                s0 = min(x[2] for x in grp) - 1
                s1 = max(x[3] for x in grp)
                at = sum(x[5] for x in grp)
                pid = sum(x[4] * x[5] for x in grp) / at
                kind = ("full_copy" if qcov >= a.full_qcov else
                        "partial_copy" if qcov >= a.part_qcov else "fragment")
                status, hocc, hbag = "new_copy", "", ""
                for (c0, c1, oid, cb) in cat.get((acc, contig), ()):
                    ov = min(s1, c1) - max(s0, c0)
                    if ov <= 0: continue
                    if ov >= a.recip * (s1 - s0) and ov >= a.recip * (c1 - c0):
                        status, hocc, hbag = "catalogued", oid, cb; break
                n[kind] += 1; n[kind + "|" + status] += 1
                out.write("\t".join(str(v) for v in [
                    q, bag, acc, species.get(acc, "?"), contig, s0, s1, strand,
                    len(grp), round(qcov, 4), round(pid, 2), kind, status, hocc, hbag,
                    "" if status == "new_copy" else ("same" if hbag == bag else "other")
                ]) + "\n")

    summ = {"hsps": nh, "copies": sum(n[k] for k in ("full_copy", "partial_copy", "fragment")),
            **{k: n[k] for k in sorted(n)}}
    with open(a.out + ".json", "w") as f: json.dump(summ, f, indent=2)
    sys.stderr.write(json.dumps(summ) + "\n")


if __name__ == "__main__":
    main()
