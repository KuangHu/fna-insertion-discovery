#!/usr/bin/env python3
"""Consolidate every insertion call into ONE folder: sequence, empty site, position.

WHY THIS EXISTS. The three things you need per insertion live in three places,
at two different grains, under two different key formats:

  position        locations/insertion_locations.tsv   (event x genome x role)
  insert sequence catalogue_v3/<sp>_inserts.fna       (>E000001|sa_k12)
  empty sequence  catalogue_v3/<sp>_targets.fna       (>E000001|sa_k12)

database_v3 holds the same sequences under the nicer key >abaum.E000001 but
EXCLUDES sa_k12 (4,661 events). catalogue_v3 is therefore the source here,
because it is the only one that covers all 109,709 events. Both key formats are
normalised to db_id = <species>.<event_id>.

THREE MODES, BECAUSE ONE IS EXPENSIVE.
  tables   joins the catalogue to the coordinates and emits the per-event table
           plus the two sequence FASTAs. No genome reads; a few minutes.
  context  slices [flank][empty junction][flank] out of each empty assembly.
           This needs 11,192 genome reads, so it groups events by genome, reads
           each one ONCE, and sharding is on the genome, never on the event --
           sharding on events would re-read the same assembly in every shard.
  merge    concatenates the context shards and checks the count.

THE CONTEXT IS CUT FROM THE EMPTY GENOME, NOT A CARRIER. coord_frame names the
pre-event assembly; junction_L_abs/junction_R_abs are positions in it. The cut
therefore shows the site as it looks WITHOUT the insertion -- upstream and
downstream abutting -- which is the thing you compare the carrier against.
"""
import argparse
import collections
import os
import sqlite3
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
OUT = W + "/insertions_all"
CAT = W + "/catalogue_v3"
DB = W + "/locations/insertions.sqlite"

EVENT_COLS = [
    "db_id", "species", "event_id", "run", "panel", "locus_id",
    "coord_frame", "coord_contig", "coord_strand",
    "junction_L_abs", "junction_R_abs", "locus_start", "locus_end",
    "inserted_len", "insert_md5", "target_len", "insertion_point_offset",
    "junction_overlap_bp", "target_bases_lost", "decomposition_method",
    "placement_status", "S1_structurally_clean",
    "insert_copies_in_genome", "insert_distinct_loci",
    "M1_within_genome", "M2_cross_event", "mobility_positive",
    "n_loci_in_event", "n_empty_genomes", "n_derived_genomes",
    "n_alleles", "n_genomes_placed",
    "n_carriers_listed", "n_empty_listed",
    "carriers_derived", "carriers_empty",
]


def species_list():
    return sorted(f[:-len("_events.tsv")] for f in os.listdir(CAT)
                  if f.endswith("_events.tsv"))


def load_locations():
    """db_id -> per-event coordinate meta, from the SQLite index."""
    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    con.row_factory = sqlite3.Row
    meta = {}
    q = con.execute(
        "SELECT db_id, run, panel, locus_id, coord_frame, coord_contig,"
        "  coord_strand, junction_L_abs, junction_R_abs, locus_start,"
        "  locus_end, n_empty_genomes, n_derived_genomes, n_alleles,"
        "  n_genomes_placed FROM loc GROUP BY db_id")
    for r in q:
        meta[r["db_id"]] = dict(r)
    return meta


def iter_fasta(path):
    key, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if key is not None:
                    yield key, hdr, "".join(buf)
                hdr = line[1:].rstrip()
                key = hdr.split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if key is not None:
        yield key, hdr, "".join(buf)


def norm_key(key, sp):
    """'E000001|sa_k12' or 'abaum.E000001' -> 'sa_k12.E000001'."""
    if "|" in key:
        ev, s = key.split("|", 1)
        return "%s.%s" % (s, ev)
    if "." in key:
        return key
    return "%s.%s" % (sp, key)


def mode_tables(a):
    os.makedirs(OUT, exist_ok=True)
    meta = load_locations()
    sys.stderr.write("locations: %d events\n" % len(meta))

    f_tsv = open(OUT + "/insertions.tsv", "w")
    f_tsv.write("\t".join(EVENT_COLS) + "\n")
    f_ins = open(OUT + "/insert.fna", "w")
    f_emp = open(OUT + "/empty_site.fna", "w")

    n_row = n_ins = n_emp = n_nocoord = 0
    for sp in species_list():
        ev_path = "%s/%s_events.tsv" % (CAT, sp)
        with open(ev_path) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(cols)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) != len(cols):
                    continue
                db_id = "%s.%s" % (sp, f[ix["event_id"]])
                m = meta.get(db_id)
                if m is None:
                    n_nocoord += 1
                    m = {}
                cd = f[ix["carriers_derived"]]
                ce = f[ix["carriers_empty"]]
                row = {
                    "db_id": db_id, "species": sp, "event_id": f[ix["event_id"]],
                    "n_carriers_listed": len([x for x in cd.split(",") if x]),
                    "n_empty_listed": len([x for x in ce.split(",") if x]),
                }
                for c in EVENT_COLS:
                    if c in row:
                        continue
                    if c in m:
                        row[c] = m[c]
                    elif c in ix:
                        row[c] = f[ix[c]]
                    else:
                        row[c] = "."
                f_tsv.write("\t".join(str(row[c]) for c in EVENT_COLS) + "\n")
                n_row += 1

        for tag, fh_out in (("inserts", f_ins), ("targets", f_emp)):
            p = "%s/%s_%s.fna" % (CAT, sp, tag)
            if not os.path.exists(p):
                continue
            for key, hdr, seq in iter_fasta(p):
                k = norm_key(key, sp)
                rest = hdr.split(None, 1)
                rest = rest[1] if len(rest) > 1 else ""
                fh_out.write(">%s %s\n" % (k, rest))
                for i in range(0, len(seq), 80):
                    fh_out.write(seq[i:i + 80] + "\n")
                if tag == "inserts":
                    n_ins += 1
                else:
                    n_emp += 1
        sys.stderr.write("  %-16s rows=%d\n" % (sp, n_row))

    for fh in (f_tsv, f_ins, f_emp):
        fh.close()
    print("insertions.tsv  %d rows" % n_row)
    print("insert.fna      %d records" % n_ins)
    print("empty_site.fna  %d records" % n_emp)
    print("events with no coordinates: %d" % n_nocoord)
    if n_row != n_ins or n_row != n_emp:
        sys.exit("FATAL: table/FASTA record counts disagree "
                 "(%d/%d/%d)" % (n_row, n_ins, n_emp))


# ------------------------------------------------------------------ context

_MAN = {}
MANIFEST = [W + "/protctx/all_fna.txt"]   # replaced by --manifest


def genome_path(acc):
    if not _MAN:
        for line in open(MANIFEST[0]):
            p = line.strip()
            if not p:
                continue
            b = os.path.basename(p)
            for e in (".fna", ".fa", ".fasta"):
                if b.endswith(e):
                    b = b[:-len(e)]
                    break
            _MAN[b] = p
    return _MAN.get(acc)


def load_genome(path):
    out, k, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if k is not None:
                    out[k] = "".join(buf)
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if k is not None:
        out[k] = "".join(buf)
    return out


def mode_context(a):
    os.makedirs(OUT + "/context_parts", exist_ok=True)
    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    by_genome = collections.defaultdict(list)
    q = con.execute(
        "SELECT db_id, coord_frame, coord_contig, coord_strand,"
        "  junction_L_abs, junction_R_abs, inserted_len, insert_md5,"
        "  n_empty_genomes, n_derived_genomes FROM loc"
        "  WHERE junction_L_abs NOT IN ('.','') GROUP BY db_id")
    for r in q:
        by_genome[r[1]].append(r)

    genomes = sorted(by_genome)
    mine = [g for i, g in enumerate(genomes) if i % a.nshard == a.shard]
    sys.stderr.write("shard %d/%d: %d genomes, %d events\n"
                     % (a.shard, a.nshard, len(mine),
                        sum(len(by_genome[g]) for g in mine)))

    out = open("%s/context_parts/part%03d.fna" % (OUT, a.shard), "w")
    n_ok = n_missing = n_short = 0
    for n, g in enumerate(mine, 1):
        p = genome_path(g)
        if not p:
            n_missing += len(by_genome[g])
            continue
        try:
            ctgs = load_genome(p)
        except (IOError, OSError):
            n_missing += len(by_genome[g])
            continue
        for (db_id, frame, contig, strand, jl, jr,
             ilen, md5, nemp, nder) in by_genome[g]:
            seq = ctgs.get(contig)
            if seq is None:
                n_missing += 1
                continue
            jl, jr = int(jl), int(jr)
            lo = max(0, jl - a.flank)
            hi = min(len(seq), jr + a.flank)
            up, site, dn = seq[lo:jl], seq[jl:jr], seq[jr:hi]
            full = (jl - lo == a.flank) and (hi - jr == a.flank)
            if not full:
                n_short += 1
            out.write(">%s frame=%s contig=%s strand=%s jL=%d jR=%d "
                      "up=%d site=%d dn=%d full=%s inserted_len=%s "
                      "insert_md5=%s n_empty=%s n_carriers=%s\n"
                      % (db_id, frame, contig, strand, jl, jr,
                         len(up), len(site), len(dn), "yes" if full else "no",
                         ilen, md5, nemp, nder))
            s = up + site + dn
            for i in range(0, len(s), 80):
                out.write(s[i:i + 80] + "\n")
            n_ok += 1
        if n % 200 == 0:
            sys.stderr.write("  %d/%d genomes, %d records\n" % (n, len(mine), n_ok))
    out.close()
    print("shard %d: %d records, %d missing, %d truncated flank"
          % (a.shard, n_ok, n_missing, n_short))


def mode_merge(a):
    parts = sorted(os.listdir(OUT + "/context_parts"))
    n = 0
    with open(OUT + "/empty_context.fna", "w") as out:
        for p in parts:
            with open(OUT + "/context_parts/" + p) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    out.write(line)
    print("empty_context.fna  %d records from %d parts" % (n, len(parts)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["tables", "context", "merge"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--manifest", help="genome list; default is the protctx corpus, "
                                       "which does NOT cover the deep_* runs")
    a = ap.parse_args()
    if a.manifest:
        MANIFEST[0] = a.manifest
    {"tables": mode_tables, "context": mode_context, "merge": mode_merge}[a.mode](a)


if __name__ == "__main__":
    main()
