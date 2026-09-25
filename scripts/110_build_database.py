#!/usr/bin/env python3
"""THE DATABASE: every (empty target site, insertion) pair, all species, one table.

This is the product the project was founded to produce. Each row is one distinct
event and carries the two sequences plus the offset that joins them:

    target_seq[:offset] + insert_seq + target_seq[offset:]  ==  the observed
    derived allele

so the empty site, the element and the exact junction are all recoverable from
the row alone. Nothing here was named, classified or matched against a library
at any point.

IDS ARE NAMESPACED BY SPECIES. `event_id` is panel-scoped within a species and
collides across them; the database key is `<species>.<event_id>`. Getting this
wrong is the eighth failure of its kind in this project, so it is done once,
here, and the row carries both parts.

WHAT A ROW ASSERTS: at this anchored site these genomes differ by an inserted
segment; `target_seq` is the empty allele verbatim; the insert goes in at
`insertion_point_offset`.

WHAT NO ROW ASSERTS: which allele is ancestral (no polarity in this pipeline);
what family the insert belongs to (no annotation); that the junction is correct
to the base (assembly tier, never read-verified); that the insert is a mobile
element (`M1`/`M2` are structural evidence, and `>=2 distinct contexts` is NOT
yet validated as insertion-specific -- see the targeted-null result).

Columns are copied through unchanged so the database cannot drift from the
per-species catalogues it is built from.

    110_build_database.py --catalogue <dir> --out database/insertions
"""
import argparse, csv, glob, os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log                                          # noqa: E402

PHYLUM = {"ecoli": "Gammaproteobacteria", "senterica": "Gammaproteobacteria",
          "kpneu": "Gammaproteobacteria", "abaum": "Gammaproteobacteria",
          "paeruginosa": "Gammaproteobacteria",
          "efaecium": "Firmicutes", "efaecalis": "Firmicutes",
          "sa_k24": "Firmicutes", "sa_k12": "Firmicutes",
          "spneumoniae": "Firmicutes", "lmonocytogenes": "Firmicutes",
          "bsubtilis": "Firmicutes",
          "mtuberculosis": "Actinobacteria",
          "ngonorrhoeae": "Betaproteobacteria", "bpertussis": "Betaproteobacteria",
          "hpylori": "Campylobacterota", "cjejuni": "Campylobacterota",
          "tpallidum": "Spirochaetes", "ctrachomatis": "Chlamydiae",
          "mpneumoniae": "Mollicutes"}
NAME = {"ecoli": "Escherichia coli", "senterica": "Salmonella enterica",
        "kpneu": "Klebsiella pneumoniae", "abaum": "Acinetobacter baumannii",
        "paeruginosa": "Pseudomonas aeruginosa", "efaecium": "Enterococcus faecium",
        "efaecalis": "Enterococcus faecalis", "sa_k24": "Staphylococcus aureus",
        "sa_k12": "Staphylococcus aureus (k12)",
        "spneumoniae": "Streptococcus pneumoniae",
        "lmonocytogenes": "Listeria monocytogenes", "bsubtilis": "Bacillus subtilis",
        "mtuberculosis": "Mycobacterium tuberculosis",
        "ngonorrhoeae": "Neisseria gonorrhoeae", "bpertussis": "Bordetella pertussis",
        "hpylori": "Helicobacter pylori", "cjejuni": "Campylobacter jejuni",
        "tpallidum": "Treponema pallidum", "ctrachomatis": "Chlamydia trachomatis",
        "mpneumoniae": "Mycoplasmoides pneumoniae"}

COLS = ["db_id", "species", "species_name", "phylum", "event_id", "locus_id",
        "target_seq", "target_len", "insertion_point_offset",
        "insert_seq", "inserted_len", "insert_md5",
        "junction_overlap_bp", "target_bases_lost", "decomposition_method",
        "placement_status", "S1_structurally_clean",
        "insert_copies_in_genome", "insert_distinct_loci",
        "M1_within_genome", "M2_cross_event", "mobility_positive",
        "n_loci_in_event", "carriers_derived", "carriers_empty"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude", nargs="*", default=["sa_k12"],
                    help="species keys to leave out of the merged table; sa_k12 "
                         "is a second arm of S. aureus for the k experiment and "
                         "would double-count it")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.catalogue, "*_events.tsv")))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n_tot = 0
    per = []
    miss = collections.Counter()
    with open(args.out + ".tsv", "w") as out, \
            open(args.out + "_targets.fna", "w") as ft, \
            open(args.out + "_inserts.fna", "w") as fi:
        out.write("\t".join(COLS) + "\n")
        for f in files:
            sp = os.path.basename(f)[:-len("_events.tsv")]
            if sp in args.exclude:
                log("skipping %s (excluded)" % sp)
                continue
            rows = list(csv.DictReader(open(f), delimiter="\t"))
            if not rows:
                continue
            n_ok = 0
            for r in rows:
                # namespaced key: event_id is panel-scoped WITHIN a species and
                # collides across them
                db = "%s.%s" % (sp, r["event_id"])
                tgt, ins = r.get("target_seq", ""), r.get("insert_seq", "")
                if not tgt or not ins:
                    miss["no_sequence"] += 1
                    continue
                try:
                    off = int(r["insertion_point_offset"])
                except (ValueError, KeyError):
                    miss["bad_offset"] += 1
                    continue
                if off < 0 or off > len(tgt):
                    miss["offset_out_of_range"] += 1
                    continue
                d = {c: r.get(c, ".") for c in COLS}
                d.update({"db_id": db, "species": sp,
                          "species_name": NAME.get(sp, sp),
                          "phylum": PHYLUM.get(sp, "?")})
                out.write("\t".join(str(d[c]) for c in COLS) + "\n")
                ft.write(">%s offset=%d len=%d overlap=%s\n%s\n"
                         % (db, off, len(tgt), r.get("junction_overlap_bp", "."), tgt))
                fi.write(">%s len=%d md5=%s\n%s\n"
                         % (db, len(ins), r.get("insert_md5", "."), ins))
                n_ok += 1
            per.append((sp, len(rows), n_ok))
            n_tot += n_ok

    log("=" * 78)
    log("INSERTION DATABASE")
    log("")
    log("  %-26s %-20s %9s %9s" % ("species", "phylum", "events", "in db"))
    log("  " + "-" * 68)
    for sp, n, ok in sorted(per, key=lambda x: -x[2]):
        log("  %-26s %-20s %9d %9d"
            % (NAME.get(sp, sp)[:26], PHYLUM.get(sp, "?")[:20], n, ok))
    log("  " + "-" * 68)
    log("  %-26s %-20s %9s %9d"
        % ("TOTAL  (%d species)" % len(per), "%d phyla"
           % len({PHYLUM.get(s, "?") for s, _, _ in per}), "", n_tot))
    if miss:
        log("")
        log("  rows dropped: %s" % dict(miss))
    log("")
    log("  Each row regenerates its derived allele:")
    log("    target_seq[:offset] + insert_seq + target_seq[offset:]")
    log("=" * 78)
    for s in (".tsv", "_targets.fna", "_inserts.fna"):
        log("wrote %s%s" % (args.out, s))
    return 0


if __name__ == "__main__":
    sys.exit(main())
