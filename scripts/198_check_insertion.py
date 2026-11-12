#!/usr/bin/env python3
"""Look up an insertion by accession+position (or id) and show it in context.

Answers, for one site: what is inserted there, what lies upstream and
downstream, where the empty version of the same site is, and how many
assemblies carry that empty version.

THE POSITION IS IN THE EMPTY GENOME. A query gives an accession and a
coordinate; that coordinate is matched against junction_L_abs/junction_R_abs,
which are positions in the PRE-EVENT assembly named by coord_frame. Querying a
carrier's accession with a coordinate read off the carrier will not match --
carrier-side coordinates were never computed. Query by --db-id, or by the
empty genome's accession. --by-carrier finds events a given accession carries,
but reports them in the empty genome's frame, which is the only frame that
exists.

THREE SEQUENCES, THREE SOURCES.
  upstream/downstream  sliced live from the empty genome's FASTA, so the window
                       is whatever the caller asks for and never a stale copy.
  the empty site       the target interval between the two junctions.
  the insertion        from database_v3/insertions_inserts.fna, which is the
                       catalogue's own canonical insert, NOT re-derived here.
                       Re-deriving it would mean redoing the LCP/LCS
                       decomposition, where lcp_bp is a SHORT-allele offset and
                       slicing the long allele with it corrupts the result.

sa_k12 IS IN THE LOCATIONS TABLE BUT NOT IN database_v3 (4,661 events, excluded
there by design). Those events resolve to coordinates and empty-site counts but
have no stored insert sequence; the script says so rather than printing nothing.
"""
import argparse
import collections
import csv
import os
import sqlite3
import sys

csv.field_size_limit(sys.maxsize)
W = "/global/scratch/users/kh36969/fna_ins_discovery"
DB = W + "/locations/insertions.sqlite"
LOC = W + "/locations/insertion_locations.tsv"
FNA = {"insert": W + "/database_v3/insertions_inserts.fna",
       "target": W + "/database_v3/insertions_targets.fna"}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def wrap(s, n=60, pad="    "):
    return "\n".join(pad + s[i:i + n] for i in range(0, len(s), n))


# ---------------------------------------------------------------- lookup

def rows_sqlite(a):
    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    con.row_factory = sqlite3.Row
    c = con.cursor()
    if a.db_id:
        q = c.execute("SELECT * FROM loc WHERE db_id=?", (a.db_id,))
    elif a.by_carrier:
        ids = [r[0] for r in c.execute(
            "SELECT DISTINCT db_id FROM loc WHERE genome=? AND role='carrier_derived'",
            (a.by_carrier,))]
        if not ids:
            return []
        q = c.execute("SELECT * FROM loc WHERE db_id IN (%s)"
                      % ",".join("?" * len(ids)), ids)
    else:
        # junction_L_abs is TEXT; compare numerically so the range is real
        q = c.execute(
            "SELECT * FROM loc WHERE db_id IN ("
            "  SELECT DISTINCT db_id FROM loc WHERE coord_frame=?"
            "    AND junction_L_abs NOT IN ('.','')"
            "    AND CAST(junction_L_abs AS INTEGER) BETWEEN ? AND ?)",
            (a.accession, a.position - a.window, a.position + a.window))
    return [dict(r) for r in q.fetchall()]


def rows_tsv(a):
    out = []
    for r in csv.DictReader(open(LOC), delimiter="\t"):
        if a.db_id:
            if r["db_id"] != a.db_id:
                continue
        elif a.by_carrier:
            pass  # filtered in a second pass below
        else:
            if r["coord_frame"] != a.accession:
                continue
            if r["junction_L_abs"] in (".", ""):
                continue
            jl, jr = int(r["junction_L_abs"]), int(r["junction_R_abs"])
            if not (jl - a.window <= a.position <= jr + a.window):
                continue
        out.append(r)
    if a.by_carrier:
        keep = set(r["db_id"] for r in out
                   if r["genome"] == a.by_carrier and r["role"] == "carrier_derived")
        out = [r for r in out if r["db_id"] in keep]
    return out


def fetch_seq(db_id, kind):
    """(header, sequence) for one db_id from the catalogue FASTA, or None."""
    if os.path.exists(DB):
        con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
        r = con.execute("SELECT header,off,nbytes FROM seqs WHERE db_id=? AND kind=?",
                        (db_id, kind)).fetchone()
        if not r:
            return None
        with open(FNA[kind], "rb") as fh:
            fh.seek(r[1])
            return r[0], fh.read(r[2]).decode().replace("\n", "")
    if not os.path.exists(FNA[kind]):
        return None
    with open(FNA[kind]) as fh:                       # fallback: linear scan
        hit, hdr, buf = False, None, []
        for line in fh:
            if line.startswith(">"):
                if hit:
                    return hdr, "".join(buf)
                hit = line[1:].split()[0] == db_id
                hdr = line[1:].rstrip()
                buf = []
            elif hit:
                buf.append(line.strip())
        return (hdr, "".join(buf)) if hit else None


# ---------------------------------------------------------------- genome

_MAN = {}


def genome_path(acc):
    if not _MAN:
        for line in open(W + "/protctx/all_fna.txt"):
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


def load_contig(path, want):
    k, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if k == want:
                    return "".join(buf)
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    return "".join(buf) if k == want else None


# ---------------------------------------------------------------- report

def report(db_id, e, a):
    m = e["meta"]
    print("=" * 72)
    print("%s   %s   locus %s" % (db_id, m["species"], m["locus_id"]))
    print("  run/panel    %s  %s" % (m["run"], m["panel"]))
    print("  inserted_len %s   insert_md5 %s" % (m["inserted_len"], m["insert_md5"]))
    print()
    print("  EMPTY SITE (every coordinate below is in this assembly)")
    print("    genome   %s" % m["coord_frame"])
    print("    contig   %s  strand %s" % (m["coord_contig"], m["coord_strand"]))
    print("    junction %s .. %s" % (m["junction_L_abs"], m["junction_R_abs"]))
    print("    locus    %s .. %s" % (m["locus_start"], m["locus_end"]))
    print()
    print("  EMPTY SITES  (assemblies showing this site with no insertion): %s"
          % m["n_empty_genomes"])
    print("  CARRIERS     (assemblies with the insertion):                  %s"
          % m["n_derived_genomes"])
    print("  alleles at this locus: %s   genomes placed: %s"
          % (m["n_alleles"], m["n_genomes_placed"]))
    print()
    lim = None if a.all_genomes else 8
    for lab, key in (("empty  ", "empty"), ("carrier", "derived")):
        g = [x["genome"] for x in e[key]]
        shown = g if lim is None else g[:lim]
        print("  %s assemblies (%d): %s%s"
              % (lab, len(g), ", ".join(shown),
                 "" if lim is None or len(g) <= lim else "  ... (--all-genomes)"))

    if not a.seq:
        return
    if m["junction_L_abs"] in (".", ""):
        print("\n  (no junction coordinates for this event)")
        return
    p = genome_path(m["coord_frame"])
    if not p:
        print("\n  (assembly %s not in the corpus manifest)" % m["coord_frame"])
        return
    ctg = load_contig(p, m["coord_contig"])
    if not ctg:
        print("\n  (contig %s not found in %s)" % (m["coord_contig"], p))
        return
    jl, jr = int(m["junction_L_abs"]), int(m["junction_R_abs"])
    up = ctg[max(0, jl - a.flank):jl]
    site = ctg[jl:jr]
    dn = ctg[jr:jr + a.flank]

    print()
    print("  --- EMPTY ALLELE, %d bp of flank from %s ---" % (a.flank, m["coord_frame"]))
    print("  UPSTREAM (%d bp)" % len(up))
    print(wrap(up))
    print("  SITE between the junctions (%d bp)%s"
          % (len(site), "   <- the insertion goes here" if not site else ""))
    if site:
        print(wrap(site))
    print("  DOWNSTREAM (%d bp)" % len(dn))
    print(wrap(dn))

    got = fetch_seq(db_id, "insert")
    print()
    if got:
        print("  --- THE INSERTION (%d bp, catalogue canonical) ---" % len(got[1]))
        print("  %s" % got[0])
        print(wrap(got[1] if a.full_insert else got[1][:600]))
        if not a.full_insert and len(got[1]) > 600:
            print("    ... %d bp more (--full-insert)" % (len(got[1]) - 600))
    else:
        print("  --- THE INSERTION ---")
        print("  not in database_v3 (sa_k12 is excluded there by design);"
              " inserted_len %s, md5 %s" % (m["inserted_len"], m["insert_md5"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", help="EMPTY genome accession")
    ap.add_argument("--position", type=int, help="coordinate on that genome")
    ap.add_argument("--window", type=int, default=2000, help="slop around --position")
    ap.add_argument("--db-id", help="look up directly, e.g. ecoli.E004606")
    ap.add_argument("--by-carrier", help="events this accession CARRIES")
    ap.add_argument("--flank", type=int, default=1000, help="bp of context to show")
    ap.add_argument("--seq", action="store_true", help="print sequences")
    ap.add_argument("--full-insert", action="store_true", help="whole insert, not 600bp")
    ap.add_argument("--all-genomes", action="store_true", help="list every assembly")
    ap.add_argument("--max-events", type=int, default=20)
    a = ap.parse_args()
    if not (a.db_id or a.by_carrier or (a.accession and a.position is not None)):
        sys.exit("give --db-id, --by-carrier, or --accession with --position")

    rows = rows_sqlite(a) if os.path.exists(DB) else rows_tsv(a)
    ev = collections.OrderedDict()
    for r in rows:
        e = ev.setdefault(r["db_id"], {"meta": r, "derived": [], "empty": []})
        (e["derived"] if r["role"] == "carrier_derived" else e["empty"]).append(r)

    if not ev:
        print("no insertion found for that query")
        print("  NOTE: --position is matched against the EMPTY genome's frame"
              " (coord_frame). A carrier accession will not match;"
              " use --by-carrier for that.")
        return

    print("%d event(s) match" % len(ev))
    for i, (db_id, e) in enumerate(ev.items()):
        if i >= a.max_events:
            print("\n... %d more (--max-events)" % (len(ev) - a.max_events))
            break
        report(db_id, e, a)


if __name__ == "__main__":
    main()
