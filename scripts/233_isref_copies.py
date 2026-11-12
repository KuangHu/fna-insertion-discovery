#!/usr/bin/env python3
"""Turn the reference-guided blastn sweep into copies, and say which copies the
discovery pipeline already knew about.

WHAT THIS IS.  ISseeker's first step is "blast the IS against a collection of
contigs"; its second is "extract the flanks and map them to a reference to
prove an insertion".  This script is step one plus the bookkeeping ISseeker
does not do, because ISseeker has no catalogue to compare against and we do.

WHAT THIS IS NOT.  A copy is not an insertion event.  An element present in
every genome at the same locus is one ancestral copy, not N insertions.  The
event test needs an empty allele and lives in the next script.  Every count
here is therefore labelled `copy`, never `event`.

Chaining.  blastn returns HSPs, and one real copy of a 1.3 kb IS can come back
as three HSPs split by SNP clusters.  HSPs on one subject are chained into a
copy when they are collinear in query and subject and the subject gap is less
than the query length -- two copies of the same IS 400 bp apart on one contig
would otherwise be merged into one, and two copies 4 kb apart would otherwise
be reported as six fragments.
"""
import argparse, collections, csv, json, os, sys

def merged_len(iv):
    iv = sorted(iv); tot = 0; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce: ce = max(ce, e)
        else: tot += ce - cs; cs, ce = s, e
    return tot + ce - cs

def chain(hsps, qlen):
    """hsps: list of (qs,qe,ss,se,pid,alen,strand). -> list of copies."""
    out = []
    for strand in ("+", "-"):
        h = sorted([x for x in hsps if x[6] == strand], key=lambda x: x[2])
        cur = []
        for x in h:
            if not cur:
                cur = [x]; continue
            gap = x[2] - max(y[3] for y in cur)
            if gap <= qlen:
                cur.append(x)
            else:
                out.append((strand, cur)); cur = [x]
        if cur: out.append((strand, cur))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--occurrences", required=True)
    ap.add_argument("--occbags", required=True, help="strict_occurrence_bags.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--full-qcov", type=float, default=0.90)
    ap.add_argument("--part-qcov", type=float, default=0.50)
    ap.add_argument("--recip", type=float, default=0.50,
                    help="reciprocal overlap to call a copy already catalogued")
    a = ap.parse_args()

    qmeta = {r["bag_id"]: r for r in csv.DictReader(open(a.queries), delimiter="\t")}
    sys.stderr.write("queries %d\n" % len(qmeta))

    species = {}
    for r in csv.DictReader(open(a.corpus), delimiter="\t"):
        species[r["assembly"]] = r["species"]
    sys.stderr.write("corpus assemblies %d\n" % len(species))

    # occurrence -> bag
    occbag = {}
    with open(a.occbags) as f:
        rd = csv.DictReader(f, delimiter="\t")
        # 228 writes the bag id into the cds_cluster_id column: bag_id IS the
        # cluster id, so accept either spelling rather than silently keying on
        # a column that does not exist.
        bcol = ("bag_id" if "bag_id" in rd.fieldnames else
                "cds_cluster_id" if "cds_cluster_id" in rd.fieldnames else None)
        ocol = "occurrence_id" if "occurrence_id" in rd.fieldnames else None
        if not (bcol and ocol):
            sys.exit("FATAL: occbags needs occurrence_id and bag_id/cds_cluster_id; got %s"
                     % rd.fieldnames)
        for r in rd: occbag[r[ocol]] = r[bcol]
    sys.stderr.write("occurrences with a bag %d\n" % len(occbag))

    # catalogue intervals, keyed by (assembly, contig)
    cat = collections.defaultdict(list)
    n_occ = 0
    for r in csv.DictReader(open(a.occurrences), delimiter="\t"):
        k = (r["carrier_assembly_accession_version"], r["carrier_fasta_record_id"])
        cat[k].append((int(r["insert_start0"]), int(r["insert_end0"]),
                       r["occurrence_id"], occbag.get(r["occurrence_id"], "")))
        n_occ += 1
    for k in cat: cat[k].sort()
    sys.stderr.write("catalogue occurrences %d over %d (assembly,contig)\n" % (n_occ, len(cat)))

    # ---- read hits, group by (query, subject)
    by = collections.defaultdict(list)
    nh = 0
    for line in open(a.hits):
        f = line.rstrip("\n").split("\t")
        q, s = f[0], f[1]
        pid = float(f[2]); alen = int(f[3])
        qs, qe = int(f[6]), int(f[7]); ss, se = int(f[8]), int(f[9])
        strand = "+" if se >= ss else "-"
        if se < ss: ss, se = se, ss
        by[(q, s)].append((qs, qe, ss, se, pid, alen, strand))
        nh += 1
    sys.stderr.write("HSPs %d over %d (query,subject) pairs\n" % (nh, len(by)))

    copies = []
    for (q, s), hsps in by.items():
        qlen = int(qmeta[q]["insert_len_nt"])
        acc, _, contig = s.partition("|")
        for strand, grp in chain(hsps, qlen):
            qcov = merged_len([(x[0], x[1]) for x in grp]) / qlen
            s0 = min(x[2] for x in grp) - 1          # blast is 1-based inclusive
            s1 = max(x[3] for x in grp)
            alen_tot = sum(x[5] for x in grp)
            pid = sum(x[4] * x[5] for x in grp) / alen_tot
            kind = ("full_copy" if qcov >= a.full_qcov else
                    "partial_copy" if qcov >= a.part_qcov else "fragment")
            # already in the catalogue?
            status, hit_occ, hit_bag = "new_copy", "", ""
            for (c0, c1, oid, bag) in cat.get((acc, contig), ()):
                ov = min(s1, c1) - max(s0, c0)
                if ov <= 0: continue
                if ov >= a.recip * (s1 - s0) and ov >= a.recip * (c1 - c0):
                    status = "catalogued"; hit_occ, hit_bag = oid, bag
                    break
            copies.append(dict(
                bag_id=q, assembly=acc, contig=contig, species=species.get(acc, "?"),
                s0=s0, s1=s1, strand=strand, n_hsp=len(grp),
                qcov=round(qcov, 4), pident=round(pid, 2), copy_kind=kind,
                catalogue_status=status, catalogue_occurrence=hit_occ,
                catalogue_bag=hit_bag,
                same_bag=("" if status == "new_copy" else
                          ("same" if hit_bag == q else "other"))))
    sys.stderr.write("copies %d\n" % len(copies))

    cols = ["bag_id","assembly","species","contig","s0","s1","strand","n_hsp",
            "qcov","pident","copy_kind","catalogue_status","catalogue_occurrence",
            "catalogue_bag","same_bag"]
    with open(a.out + "_copies.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for c in sorted(copies, key=lambda c: (c["bag_id"], c["assembly"], c["contig"], c["s0"])):
            out.write("\t".join(str(c[k]) for k in cols) + "\n")

    # ---- per-bag summary
    agg = collections.defaultdict(lambda: collections.Counter())
    gen = collections.defaultdict(set); spp = collections.defaultdict(set)
    for c in copies:
        b = c["bag_id"]; agg[b][c["copy_kind"]] += 1
        agg[b]["all"] += 1
        if c["copy_kind"] == "full_copy":
            agg[b]["full_" + c["catalogue_status"]] += 1
            if c["catalogue_status"] == "catalogued":
                agg[b]["full_cat_" + c["same_bag"]] += 1
            gen[b].add(c["assembly"]); spp[b].add(c["species"])
    scols = ["bag_id","n_independent_sites","band","insert_len_nt","catalytic_families",
             "copies_all","full_copies","partial_copies","fragments",
             "full_catalogued","full_cat_same_bag","full_cat_other_bag","full_new",
             "genomes_with_full_copy","species_with_full_copy"]
    with open(a.out + "_by_bag.tsv", "w") as out:
        out.write("\t".join(scols) + "\n")
        for b, r in sorted(qmeta.items()):
            A = agg[b]
            out.write("\t".join(str(x) for x in [
                b, r["n_independent_sites"], r["band"], r["insert_len_nt"],
                r["catalytic_families"], A["all"], A["full_copy"], A["partial_copy"],
                A["fragment"], A["full_catalogued"], A["full_cat_same"],
                A["full_cat_other"], A["full_new_copy"],
                len(gen[b]), len(spp[b])]) + "\n")

    tot = collections.Counter()
    for c in copies:
        tot[c["copy_kind"]] += 1
        if c["copy_kind"] == "full_copy":
            tot["full_" + c["catalogue_status"]] += 1
            if c["catalogue_status"] == "catalogued": tot["full_cat_" + c["same_bag"]] += 1
    print(json.dumps({"hsps": nh, "copies": len(copies), **tot}, indent=2))

main()
