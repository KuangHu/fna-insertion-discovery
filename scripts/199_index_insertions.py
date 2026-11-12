#!/usr/bin/env python3
"""Index the insertion-locations table into SQLite for instant lookup.

WHY. insertion_locations.tsv is 242MB / 979k rows. A linear scan costs ~5s per
query, which is fine once and unusable for a tool meant to be run repeatedly.
SQLite with two indexes answers the same question in milliseconds and needs no
server.

TWO QUERY PATHS, SO TWO INDEXES.
  by id        -- db_id, e.g. ecoli.E004606
  by position  -- (coord_frame, junction_L_abs). A positional query asks for
                  everything near a coordinate, so the index is on the genome
                  plus the left junction and the range scan happens inside
                  SQLite rather than in Python.

THE INSERT SEQUENCES ARE INDEXED BY OFFSET, NOT COPIED. insertions_inserts.fna
and insertions_targets.fna are 166MB and 57MB. Storing a byte offset per record
and seeking on demand keeps the index small and avoids a second copy of the
sequence data going stale against the FASTA it came from.
"""
import os
import sqlite3
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
LOC = W + "/locations/insertion_locations.tsv"
DB = W + "/locations/insertions.sqlite"
FNA = {"insert": W + "/database_v3/insertions_inserts.fna",
       "target": W + "/database_v3/insertions_targets.fna"}


def index_fasta(path):
    """db_id -> (header, seq_offset, seq_nbytes). One pass, binary, no parsing."""
    out = {}
    with open(path, "rb") as fh:
        key, hdr, start = None, None, None
        pos = 0
        for line in fh:
            n = len(line)
            if line.startswith(b">"):
                if key is not None:
                    out[key] = (hdr, start, pos - start)
                hdr = line[1:].decode().rstrip()
                key = hdr.split()[0]
                start = pos + n
            pos += n
        if key is not None:
            out[key] = (hdr, start, pos - start)
    return out


def main():
    if not os.path.exists(LOC):
        sys.exit("missing %s" % LOC)
    if os.path.exists(DB):
        os.unlink(DB)
    con = sqlite3.connect(DB)
    c = con.cursor()
    c.execute("PRAGMA journal_mode=OFF")
    c.execute("PRAGMA synchronous=OFF")

    with open(LOC) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        c.execute("CREATE TABLE loc (%s)" % ",".join('"%s" TEXT' % x for x in cols))
        ins = "INSERT INTO loc VALUES (%s)" % ",".join("?" * len(cols))
        buf, n = [], 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(cols):
                continue
            buf.append(f)
            n += 1
            if len(buf) >= 50000:
                c.executemany(ins, buf)
                buf = []
        if buf:
            c.executemany(ins, buf)
    print("loc rows %d" % n)

    # byte offsets into the two FASTAs, so sequences are seeked not stored
    c.execute("CREATE TABLE seqs (db_id TEXT, kind TEXT, header TEXT, "
              "off INTEGER, nbytes INTEGER)")
    for kind, path in FNA.items():
        if not os.path.exists(path):
            print("skip %s (absent)" % path)
            continue
        idx = index_fasta(path)
        c.executemany("INSERT INTO seqs VALUES (?,?,?,?,?)",
                      [(k, kind, v[0], v[1], v[2]) for k, v in idx.items()])
        print("%-7s %d records" % (kind, len(idx)))

    print("building indexes ...")
    c.execute("CREATE INDEX i_dbid ON loc(db_id)")
    c.execute("CREATE INDEX i_pos  ON loc(coord_frame, junction_L_abs)")
    c.execute("CREATE INDEX i_gen  ON loc(genome)")
    c.execute("CREATE INDEX i_seq  ON seqs(db_id, kind)")
    con.commit()
    con.close()
    print("wrote %s  (%.0f MB)" % (DB, os.path.getsize(DB) / 1e6))


if __name__ == "__main__":
    main()
