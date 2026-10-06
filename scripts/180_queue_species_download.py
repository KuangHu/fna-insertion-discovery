#!/usr/bin/env python3
"""Build a download queue for one species: usable assemblies only.

USABLE = assembly_level Complete Genome or Chromosome. The census ranks on
this and never on total, because Arm B compares whole assemblies and a
contig-level pool contributes fragmentation rather than alleles. Salmonella is
the standing example: largest species in GenBank at 626,212 assemblies, 0.57%
usable.

MATCHED ON species_taxid, NOT organism_name. organism_name carries strain
suffixes ("Escherichia coli O157:H7 str. Sakai"), so matching on it shatters a
species across hundreds of distinct strings and silently loses most of the
pool. taxid is the stable key.

EXCLUDES anomalous and suppressed records -- NCBI flags assemblies withdrawn
for contamination or misidentification, and those would enter the panel as
ordinary genomes.

Writes <out>/<species>.tsv as accession<TAB>url, which is what
tools/download_on_transfer_node.sh consumes. Downloads then run BY HAND on the
transfer node: cf1 is OverSubscribe=EXCLUSIVE so a download submitted there is
billed a whole 64-core node while using no CPU -- measured at 225 CPU-hours
for K. pneumoniae alone.
"""
import argparse
import collections
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--taxid", required=True, type=int,
                    help="species_taxid, not taxid")
    ap.add_argument("--name", required=True, help="short label for the queue file")
    ap.add_argument("--out", required=True)
    ap.add_argument("--levels", default="Complete Genome,Chromosome")
    a = ap.parse_args()

    levels = {x.strip() for x in a.levels.split(",")}
    rows, st = [], collections.Counter()
    with open(a.summary) as fh:
        hdr = None
        for line in fh:
            if line.startswith("#assembly_accession"):
                hdr = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not hdr:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            r = dict(zip(hdr, f))
            if r.get("species_taxid") != str(a.taxid):
                continue
            st["taxid_match"] += 1
            if r.get("assembly_level") not in levels:
                st["wrong_level"] += 1
                continue
            # withdrawn for contamination or misidentification
            if r.get("excluded_from_refseq", "").strip() not in ("", "na"):
                st["excluded_from_refseq"] += 1
                continue
            url = r.get("ftp_path", "")
            if not url or url == "na":
                st["no_ftp_path"] += 1
                continue
            acc = r["assembly_accession"]
            rows.append((acc, "%s/%s_genomic.fna.gz" % (url, os.path.basename(url))))
            st["QUEUED"] += 1

    os.makedirs(a.out, exist_ok=True)
    p = os.path.join(a.out, "%s.tsv" % a.name)
    with open(p, "w") as fh:
        for acc, url in sorted(rows):
            fh.write("%s\t%s\n" % (acc, url))
    print("%-22s taxid %-8d queued %5d   (matched %d, wrong_level %d, "
          "excluded %d, no_url %d)"
          % (a.name, a.taxid, st["QUEUED"], st["taxid_match"],
             st["wrong_level"], st["excluded_from_refseq"], st["no_ftp_path"]))
    if st["QUEUED"] == 0:
        print("  FATAL: nothing queued -- check the species_taxid")
        return 2
    print("  wrote %s" % p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
