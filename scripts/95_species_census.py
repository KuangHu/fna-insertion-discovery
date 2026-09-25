#!/usr/bin/env python3
"""Per-species census of the discovery substrate. Read-only, metadata only.

WHY THIS EXISTS: every threshold in this project was tuned on E. coli, and §9 of
NEXT_STEPS.md forbids hard-coding them elsewhere. Before any parameter can be
recalibrated for a second species, the substrate has to be counted: how many
genomes exist per species, and -- the number that actually matters -- how many
are at an assembly level Arm B can use.

ARM B NEEDS COMPLETE OR CHROMOSOME. The graph-bubble arm compares whole
assemblies; a 385,370-contig-level pool contributes fragmentation, not alleles.
E. coli has 423,319 GenBank entries but only 8,958 complete+chromosome, so
ranking species by total genome count would pick the wrong pilot. Species are
therefore ranked by USABLE count, with the total kept alongside so the gap is
visible.

The species key is species_taxid, not the organism_name string: organism_name
carries strain suffixes and inconsistent nomenclature, and would shatter one
species across dozens of rows.

    95_species_census.py --summary assembly_summary_genbank_bacteria.txt \\
        --out census/genbank
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

LEVELS = ("Complete Genome", "Chromosome", "Scaffold", "Contig")
USABLE = ("Complete Genome", "Chromosome")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-usable", type=int, default=50,
                    help="species below this many complete+chromosome genomes "
                         "cannot support a panel and are counted but not ranked")
    ap.add_argument("--top", type=int, default=40)
    args = ap.parse_args()

    # species_taxid -> level -> n ; and the most common organism name seen
    cnt = collections.defaultdict(collections.Counter)
    names = collections.defaultdict(collections.Counter)
    latest = collections.Counter()
    total_rows = 0
    for path in args.summary:
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 14:
                    continue
                total_rows += 1
                # version_status: skip superseded/suppressed assemblies
                if f[10] != "latest":
                    continue
                latest["latest"] += 1
                st, org, lvl = f[6], f[7], f[11]
                cnt[st][lvl] += 1
                # genus + species only: drop strain suffix
                names[st][" ".join(org.split()[:2])] += 1
    log("%d rows, %d latest, %d distinct species_taxid"
        % (total_rows, latest["latest"], len(cnt)))

    rows = []
    for st, c in cnt.items():
        usable = sum(c[l] for l in USABLE)
        rows.append({
            "species_taxid": st,
            "species": names[st].most_common(1)[0][0] if names[st] else "?",
            "complete": c["Complete Genome"], "chromosome": c["Chromosome"],
            "scaffold": c["Scaffold"], "contig": c["Contig"],
            "usable": usable, "total": sum(c.values()),
            "usable_fraction": usable / sum(c.values()) if sum(c.values()) else 0.0})
    rows.sort(key=lambda r: -r["usable"])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    cols = ["species_taxid", "species", "complete", "chromosome", "scaffold",
            "contig", "usable", "total", "usable_fraction"]
    with open(args.out + "_species_census.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.4f" % r[c]) if isinstance(r[c], float)
                               else str(r[c]) for c in cols) + "\n")

    elig = [r for r in rows if r["usable"] >= args.min_usable]
    log("=" * 96)
    log("PER-SPECIES CENSUS  (ranked by USABLE = complete + chromosome)")
    log("")
    log("  %-38s %9s %11s %9s %10s %8s"
        % ("species", "usable", "complete", "chrom", "total", "usable%"))
    log("  " + "-" * 90)
    for r in rows[:args.top]:
        log("  %-38s %9d %11d %9d %10d %7.2f%%"
            % (r["species"][:38], r["usable"], r["complete"], r["chromosome"],
               r["total"], 100 * r["usable_fraction"]))
    log("")
    log("  species with >= %d usable genomes: %d" % (args.min_usable, len(elig)))
    log("  total usable genomes across all species: %d"
        % sum(r["usable"] for r in rows))
    log("")
    log("  Ranking by TOTAL instead would pick a different pilot -- the top 5 by")
    log("  total are:")
    for r in sorted(rows, key=lambda x: -x["total"])[:5]:
        log("    %-36s total %8d  but usable %6d (%.2f%%)"
            % (r["species"][:36], r["total"], r["usable"],
               100 * r["usable_fraction"]))
    log("=" * 96)
    log("wrote %s_species_census.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
