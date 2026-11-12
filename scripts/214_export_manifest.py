#!/usr/bin/env python3
"""Stage 1 of the v1 record-format export: source_manifest.jsonl.

One line per SOURCE CONTIG that the dataset actually cites, carrying enough to
re-open the exact bytes later: assembly accession with its version, the FNA
file and its byte hash, the exact FASTA record token, the contig length, and
the hash of the normalised contig string.

WHY A CONTIG TOKEN IS NOT AN IDENTIFIER ON ITS OWN. Two different assemblies
routinely use the same token (`contig1`, `NODE_1`, ...). Resolution is always
(assembly -> file -> record), never a search for the token. The source_id
encodes that pair so nothing downstream can look a contig up by name alone.

A duplicate token INSIDE one file is a hard failure for that file, not a
"take the first" situation: it makes every citation of that token ambiguous.
"""
import argparse
import collections
import hashlib
import json
import os
import sys
import zlib

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"
V2 = W + "/inferred_empty_v2"


def source_id(acc, rec):
    return "SRC|%s|%s" % (acc, rec)


def file_sha256(p, buf=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def parse_fna(p):
    """token -> normalised sequence. Raises on a duplicate token."""
    out, k, b = {}, None, []
    with open(p) as fh:
        for l in fh:
            if l[0] == ">":
                if k is not None:
                    if k in out:
                        raise ValueError("duplicate FASTA record id: " + k)
                    out[k] = "".join(b).upper()
                k = l[1:].split()[0]
                b = []
            else:
                b.append(l.strip())
    if k is not None:
        if k in out:
            raise ValueError("duplicate FASTA record id: " + k)
        out[k] = "".join(b).upper()
    return out


def needed():
    """accession -> set of contig tokens cited by any input record."""
    need = collections.defaultdict(set)
    with open(V2 + "/occurrences.tsv") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            need[f[ix["carrier_accession"]]].add(f[ix["carrier_contig"]])
            a = f[ix["empty_accession"]]
            if a and a != ".":
                need[a].add(f[ix["empty_contig"]])
    return need


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    man = {}
    with open(SRC + "/genome_map.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                man[f[0]] = f[1]

    need = needed()
    accs = sorted(k for k in need
                  if zlib.crc32(k.encode()) % a.nshard == a.shard % a.nshard)

    out = open("%s/source_manifest_%03d.jsonl" % (a.out, a.shard), "w")
    bad = open("%s/source_failures_%03d.tsv" % (a.out, a.shard), "w")
    n = collections.Counter()
    for i, acc in enumerate(accs, 1):
        p = man.get(acc)
        if not p or not os.path.exists(p):
            bad.write("%s\t.\tgenome_not_in_manifest\n" % acc)
            n["genome_not_in_manifest"] += 1
            continue
        try:
            ctgs = parse_fna(p)
        except ValueError as e:
            bad.write("%s\t.\t%s\n" % (acc, e))
            n["duplicate_record_id_in_file"] += 1
            continue
        except (IOError, OSError) as e:
            bad.write("%s\t.\tunreadable:%s\n" % (acc, e))
            n["unreadable"] += 1
            continue
        fh256 = file_sha256(p)
        rel = os.path.relpath(p, W)
        for rec in sorted(need[acc]):
            s = ctgs.get(rec)
            if s is None:
                bad.write("%s\t%s\tcontig_not_in_file\n" % (acc, rec))
                n["contig_not_in_file"] += 1
                continue
            out.write(json.dumps({
                "source_id": source_id(acc, rec),
                # The pipeline's accessions are assembly accessions carrying
                # their .version; the per-contig public accession is a
                # different field and we do not have it, so it stays null
                # rather than being filled with the assembly's.
                "assembly_accession_version": acc,
                "assembly_local_id": None,
                "sequence_accession_version": None,
                "fna_path": rel,
                "fasta_record_id": rec,
                "contig_length_bp": len(s),
                "contig_seq_sha256": hashlib.sha256(s.encode()).hexdigest(),
                "fna_file_sha256": fh256,
                "topology": "unknown",
            }, separators=(",", ":")) + "\n")
            n["sources"] += 1
        if i % 200 == 0:
            sys.stderr.write("  %d/%d assemblies\n" % (i, len(accs)))
    out.close()
    bad.close()
    print("shard %d: assemblies %d" % (a.shard, len(accs)))
    for k, v in sorted(n.items()):
        print("   %-34s %7d" % (k, v))


if __name__ == "__main__":
    main()
