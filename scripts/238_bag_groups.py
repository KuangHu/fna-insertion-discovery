#!/usr/bin/env python3
"""Split bags into a QUESTIONABLE group and the rest.

A bag with very little evidence behind it should not be reported next to one
built from hundreds of observations, so the thin ones are separated BEFORE any
claim is made about them rather than filtered out afterwards.

TWO DIFFERENT THINNESSES, AND THEY ARE NOT THE SAME BAG SET.

  n_distinct_sequences  how many distinct insertion DNA sequences the bag has.
                        Thin here means the bag's sequence diversity is barely
                        sampled -- there is little or nothing to compare.

  n_independent_sites   how many distinct empty-site loci the bag was seen at.
                        Thin here means no evidence of MOBILITY: an element
                        seen at one locus has never been observed to move, and
                        a single structural difference can be a deletion
                        polymorphism or an assembly artifact just as easily as
                        an insertion.

A bag can be thin on one axis and thick on the other.  One distinct sequence
present at 500 sites is a perfectly ordinary IS element; 500 distinct
sequences at one site is not an element at all.  So both are computed, both
are written, and the column that decides the group is named explicitly rather
than left implicit.

This is a LABEL, not a deletion.  Questionable bags keep their queries and are
searched exactly like the rest; the group travels with the results so a reader
can see which tier a row came from.  Nothing here says a questionable bag is
wrong -- only that it has not got enough behind it to stand on its own.
"""
import argparse, csv, json, collections


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bags", required=True, help="strict_bags.tsv")
    ap.add_argument("--occbags", required=True, help="strict_occurrence_bags.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--criterion", choices=["sequences", "sites"], default="sequences",
                    help="which thinness defines the questionable group")
    ap.add_argument("--max-thin", type=int, default=2,
                    help="bags at or below this count are questionable")
    a = ap.parse_args()

    seqs = collections.defaultdict(set)
    with open(a.occbags) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["cds_cluster_id"]:
                seqs[r["cds_cluster_id"]].add(r["seq_key"])

    rows = []
    with open(a.bags) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b = r["bag_id"]
            ns = len(seqs.get(b, ()))
            nsite = int(r["n_independent_sites"])
            nocc = int(r["n_occurrences"])
            key = ns if a.criterion == "sequences" else nsite
            rows.append(dict(bag_id=b, n_distinct_sequences=ns, n_occurrences=nocc,
                             n_independent_sites=nsite, n_species=int(r["n_species"]),
                             group="questionable" if key <= a.max_thin else "standard"))

    with open(a.out, "w") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "bag_id", "group", "n_distinct_sequences", "n_occurrences",
            "n_independent_sites", "n_species"])
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["bag_id"]): w.writerow(r)

    q = [r for r in rows if r["group"] == "questionable"]
    s = [r for r in rows if r["group"] == "standard"]
    # how different the other criterion would have been, stated not hidden
    alt = "sites" if a.criterion == "sequences" else "sequences"
    altkey = "n_independent_sites" if alt == "sites" else "n_distinct_sequences"
    altq = {r["bag_id"] for r in rows if r[altkey] <= a.max_thin}
    thisq = {r["bag_id"] for r in q}
    summ = {
        "criterion": a.criterion, "max_thin": a.max_thin,
        "bags": len(rows), "questionable": len(q), "standard": len(s),
        "questionable_occurrences": sum(r["n_occurrences"] for r in q),
        "standard_occurrences": sum(r["n_occurrences"] for r in s),
        "questionable_sites": sum(r["n_independent_sites"] for r in q),
        "standard_sites": sum(r["n_independent_sites"] for r in s),
        "if_criterion_were_" + alt: len(altq),
        "bags_questionable_under_both": len(thisq & altq),
        "bags_questionable_only_under_" + a.criterion: len(thisq - altq),
        "bags_questionable_only_under_" + alt: len(altq - thisq),
        "note": "a label on the bag, not a filter; questionable bags keep their queries",
    }
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
