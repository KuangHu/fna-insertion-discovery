#!/usr/bin/env python3
"""One row per (insertion event, genome accession, role), with coordinates.

WHAT THE CATALOGUE DOES AND DOES NOT HOLD. catalogue_v3 identifies an event by
SEQUENCE -- target_seq plus insertion_point_offset -- and lists its carrier
accessions, but has no coordinate column at all. That was deliberate: events
are deduplicated on canonical sequence keys so one insertion seen in several
assemblies collapses to a single record, and an assembly-specific coordinate
cannot survive that merge. The coordinates live in the per-panel Layer 2
output, which this joins back.

WHICH GENOME THE COORDINATES DESCRIBE. junction_L_abs/junction_R_abs are
positions in the PRE-EVENT (empty) genome -- the assembly that lacks the
insertion. They locate the empty site to the base. A carrier's own coordinates
are NOT computed here: locating the insertion inside each of N carriers needs
a mapping pass per genome. Every row therefore carries coord_frame naming the
assembly its numbers refer to, so a coordinate is never read against the wrong
genome.

n_empty_genomes is the count the user asked for: how many assemblies show this
site WITHOUT the insertion.
"""
import collections
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)
W = "/global/scratch/users/kh36969/fna_ins_discovery"
OUT = sys.argv[1] if len(sys.argv) > 1 else W + "/locations/insertion_locations.tsv"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

COLS = ["db_id", "species", "event_id", "run", "panel", "locus_id",
        "genome", "role", "allele_id", "allele_len", "allele_md5",
        "coord_frame", "coord_contig", "coord_strand",
        "junction_L_abs", "junction_R_abs", "locus_start", "locus_end",
        "insert_md5", "inserted_len",
        "n_empty_genomes", "n_derived_genomes", "n_alleles",
        "n_genomes_placed", "coord_status"]

# ---- events grouped by the panel whose Layer 2 output holds their coordinates
by_panel = collections.defaultdict(list)
n_ev = 0
for f in sorted(glob.glob(W + "/catalogue_v3/*_events.tsv")):
    sp = os.path.basename(f).replace("_events.tsv", "")
    for r in csv.DictReader(open(f), delimiter="\t"):
        n_ev += 1
        by_panel[(r["run"], r["panel"])].append((sp, r))
print("events %d across %d (run, panel) groups" % (n_ev, len(by_panel)))

out = open(OUT, "w")
out.write("\t".join(COLS) + "\n")
n_row = n_loc = n_nocoord = 0
missing_panel = 0
for i, ((run, panel), evs) in enumerate(sorted(by_panel.items())):
    pb = os.path.basename(panel)
    ldir = os.path.join(W, run, "l2", pb)
    loci, alle = {}, collections.defaultdict(list)
    lp = os.path.join(ldir, "loci.tsv")
    if os.path.exists(lp):
        for r in csv.DictReader(open(lp), delimiter="\t"):
            loci[r["locus_id"]] = r
        ap = os.path.join(ldir, "alleles.tsv")
        if os.path.exists(ap):
            for r in csv.DictReader(open(ap), delimiter="\t"):
                alle[r["locus_id"]].append(r)
    else:
        missing_panel += 1

    for sp, r in evs:
        # locus_id in the catalogue is qualified "<run>/<panel>|<locus>"; the
        # Layer 2 table keys on the bare locus id
        lid = r["locus_id"].split("|")[-1]
        L = loci.get(lid)
        if L:
            n_loc += 1
        else:
            n_nocoord += 1
        g2a = {}
        for a in alle.get(lid, []):
            for g in (a.get("genomes") or "").split(","):
                g = g.strip()
                if g:
                    g2a[g] = a
        der = [x for x in (r.get("carriers_derived") or "").split(",") if x]
        emp = [x for x in (r.get("carriers_empty") or "").split(",") if x]
        base = dict(db_id="%s.%s" % (sp, r["event_id"]), species=sp,
                    event_id=r["event_id"], run=run, panel=panel,
                    locus_id=lid,
                    insert_md5=r.get("insert_md5", "."),
                    inserted_len=r.get("inserted_len", "."),
                    n_empty_genomes=len(emp), n_derived_genomes=len(der),
                    n_alleles=(L or {}).get("n_alleles", "."),
                    n_genomes_placed=(L or {}).get("n_genomes_placed", "."),
                    coord_frame=(L or {}).get("pre_event_genome", "."),
                    coord_contig=(L or {}).get("pre_event_contig", "."),
                    coord_strand=(L or {}).get("pre_event_strand", "."),
                    junction_L_abs=(L or {}).get("junction_L_abs", "."),
                    junction_R_abs=(L or {}).get("junction_R_abs", "."),
                    locus_start=(L or {}).get("start", "."),
                    locus_end=(L or {}).get("end", "."),
                    coord_status="located" if L else "no_layer2_locus")
        for g, role in [(x, "carrier_derived") for x in der] + \
                       [(x, "carrier_empty") for x in emp]:
            a = g2a.get(g)
            row = dict(base, genome=g, role=role,
                       allele_id=(a or {}).get("allele_id", "."),
                       allele_len=(a or {}).get("length", "."),
                       allele_md5=(a or {}).get("md5", "."))
            out.write("\t".join(str(row[c]) for c in COLS) + "\n")
            n_row += 1
    if (i + 1) % 500 == 0:
        print("  %d/%d panels, %d rows" % (i + 1, len(by_panel), n_row),
              file=sys.stderr)
out.close()
print("rows %d | events with coordinates %d | without %d | panels missing l2 %d"
      % (n_row, n_loc, n_nocoord, missing_panel))
print("wrote %s" % OUT)
