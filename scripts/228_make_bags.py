#!/usr/bin/env python3
"""Expand dedup-level clusters back to occurrences and build the bag table.

Clustering ran on DISTINCT insertion DNA, so a cluster is a set of sequences,
not a set of observations.  A bag is only interesting as observations, and the
expansion back out is where counts are easy to inflate.  Three different
counts are reported and never conflated:

  n_distinct_insertions  how many distinct DNA sequences carry this bag
  n_occurrences          how many (event, carrier, contig, interval) rows
  n_independent_sites    occurrences collapsed to one per EMPTY-SITE LOCUS
                         (reference contig + reference junction)

The third exists because of a mistake already made and corrected once in this
project: two alleles sitting at one locus were being counted as two sites,
which inflated every bag of size >= 5 by about 28%.  A bag's size is the
number of places it has been seen, and two alleles at one place is one place.
`n_occurrences / n_independent_sites` is printed per bag so the degree of
collapse is visible rather than buried.

SPECIES BREADTH IS A SEPARATE AXIS from size.  A bag of 400 occurrences in one
species is a different object from a bag of 400 across twelve, so both
n_species and the per-species counts are kept.

PFAM, WHERE SUPPLIED, DESCRIBES A BAG.  It never defines one.  The families
listed per bag are those found on the bag members' own proteins, with the
fraction of members carrying each, so a bag that is 3% annotated cannot be
read as a named family.
"""
import argparse, collections, json, os, statistics


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cds", required=True, help="insertion_cds_<tag>.tsv")
    ap.add_argument("--members", required=True, help="dedup_members.tsv")
    ap.add_argument("--browse", required=True)
    ap.add_argument("--pfam", default="", help="optional *_hits_on_dna.tsv or *_pfam.tsv")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()

    cds = {r["seq_key"]: r for r in read_tsv(a.cds)}
    bag_of = {k: r["cds_cluster_id"] for k, r in cds.items() if r["cds_cluster_id"]}

    members = collections.defaultdict(list)
    for r in read_tsv(a.members):
        if r["seq_key"] in bag_of:
            members[r["seq_key"]].append(r)

    need = {m["occurrence_id"] for v in members.values() for m in v}
    loci = {}
    for r in read_tsv(a.browse):
        if r["occurrence_id"] in need:
            loci[r["occurrence_id"]] = (
                r["reference_sequence_accession_version"] or r["reference_fasta_record_id"],
                r["reference_junction0"], r["carrier_assembly_accession_version"],
                r["carrier_fasta_record_id"], r["insert_start0"], r["insert_end0"])

    fam = collections.defaultdict(collections.Counter)
    if a.pfam and os.path.exists(a.pfam):
        seen = collections.defaultdict(set)
        for h in read_tsv(a.pfam):
            k = h.get("seq_key", "")
            if k in bag_of:
                seen[k].add((h["pfam_acc"].split(".")[0], h["pfam_name"]))
        for k, s in seen.items():
            for acc, nm in s:
                fam[bag_of[k]][(acc, nm)] += 1

    occ_out = open(a.out + "_occurrence_bags.tsv", "w")
    occ_out.write("occurrence_id\tevent_id\tspecies\tseq_key\torient_to_rep\t"
                  "cds_cluster_id\tdominant_cds_source\tinsert_len_nt\n")
    bags = collections.defaultdict(lambda: dict(
        keys=set(), occ=0, events=set(), sites=set(), species=collections.Counter(),
        lens=[], domaa=[], srcs=collections.Counter()))
    for k, ms in members.items():
        b = bag_of[k]
        B = bags[b]
        B["keys"].add(k)
        B["lens"].append(int(cds[k]["insert_len_nt"]))
        B["domaa"].append(int(cds[k]["dominant_cds_len_aa"]))
        B["srcs"][cds[k]["dominant_cds_source"]] += 1
        for m in ms:
            B["occ"] += 1
            B["events"].add(m["event_id"] + "|" + m["species"])
            B["species"][m["species"]] += 1
            lo = loci.get(m["occurrence_id"])
            if lo:
                B["sites"].add((lo[0], lo[1]))
            occ_out.write("\t".join([m["occurrence_id"], m["event_id"], m["species"], k,
                                     m["orient_to_rep"], b,
                                     cds[k]["dominant_cds_source"],
                                     m["length"]]) + "\n")
    occ_out.close()

    with open(a.out + "_bags.tsv", "w") as f:
        f.write("bag_id\tn_distinct_insertions\tn_occurrences\tn_events\t"
                "n_independent_sites\tocc_per_site\tn_species\ttop_species\t"
                "median_insert_nt\tmin_insert_nt\tmax_insert_nt\tmedian_dominant_aa\t"
                "dominant_cds_sources\tpfam_families\n")
        for b, B in sorted(bags.items(), key=lambda kv: -kv[1]["occ"]):
            ns = len(B["sites"]) or len(B["events"])
            fams = ";".join(f"{acc}:{nm}:{n}" for (acc, nm), n in
                            fam.get(b, collections.Counter()).most_common(5))
            f.write("\t".join(str(x) for x in [
                b, len(B["keys"]), B["occ"], len(B["events"]), len(B["sites"]),
                f"{B['occ']/max(1,ns):.2f}", len(B["species"]),
                ",".join(f"{s}:{n}" for s, n in B["species"].most_common(3)),
                int(statistics.median(B["lens"])), min(B["lens"]), max(B["lens"]),
                int(statistics.median(B["domaa"])),
                ",".join(f"{s}:{n}" for s, n in B["srcs"].most_common()),
                fams]) + "\n")

    sizes = sorted((len(B["sites"]) or len(B["events"])) for B in bags.values())
    summ = dict(bags=len(bags), insertions=len(members),
                occurrences=sum(B["occ"] for B in bags.values()),
                independent_sites=sum(len(B["sites"]) for B in bags.values()),
                bags_ge5_sites=sum(1 for s in sizes if s >= 5),
                bags_ge20_sites=sum(1 for s in sizes if s >= 20),
                singleton_bags=sum(1 for s in sizes if s == 1),
                median_bag_sites=sizes[len(sizes) // 2] if sizes else 0,
                multi_species_bags=sum(1 for B in bags.values() if len(B["species"]) > 1))
    with open(a.out + "_bags.json", "w") as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
