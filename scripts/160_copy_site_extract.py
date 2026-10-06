#!/usr/bin/env python3
"""Turn detected COPIES into bag sites by excising the element from the carrier.

Today one event yields one site. The 783,577 copies at distinct loci that
copy-finding already detects yield nothing, because a bag site needs the
PRE-INSERTION sequence and a homology hit establishes only presence.

This excises: map a known insert into genomes, and where it lands cleanly,
remove it and rejoin the flanks. The result is that carrier's pre-insertion
site at that locus.

MEASURED ACCURACY, and the exact population it was measured on (152_):

    flank-120 agreement with an OBSERVED empty allele    84.7%
      target_bases_lost == 0                             85.4%
      junction overlap 4-15                              89.3%
      junction overlap 0                                 44.3%
    when wrong: median 4 bases differ, p90 30, max 93

Two things that number is NOT:

  * It is measured on loci that HAVE an observed empty allele -- sites
    demonstrably once empty. The copies this script is for have none, may be
    ancestral, and their true accuracy is UNMEASURED and probably lower.
  * Ground truth there was a different genome's allele, so part of the 15%
    disagreement is real strain divergence rather than reconstruction error.
    84.7% is therefore a floor on accuracy and a ceiling on what can be
    measured at all. A clean test needs one genome with and without the
    element, which cannot exist.

THE GATE, and its honest status. The validated gate
(target_bases_lost == 0 AND overlap >= 1) is NOT computable for a copy -- both
fields come from a decomposition that needs the empty allele. What IS visible
from the carrier alone is a direct repeat flanking the element, so that is the
gate used here: --min-tsd. Its 100% precision in testing was confounded (it
never co-occurs with the tolerant path, so it proxied a tautology), which means
**this gate's accuracy on copies is unvalidated**. It is applied because
excision is only well defined where the boundary is, not because it is known
to work.

Every row is written with empty_site_source=tsd_reconstructed so nothing here
can be mistaken for an observed site, and a consumer can drop them wholesale.
"""
import argparse
import collections
import csv
import os
import subprocess
import sys
import tempfile

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_insert as LI
from lib.util import log                                          # noqa: E402

HALF = 60


def read_fa(p):
    out, name, buf = {}, None, []
    for line in open(p):
        if line[0] == ">":
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def tsd_len(seq, a, b, cap=40):
    """Direct repeat flanking [a,b). Carrier-side only -- no empty allele."""
    n = 0
    while n < cap and a - 1 - n >= 0 and b - 1 - n >= 0 \
            and seq[a - 1 - n] == seq[b - 1 - n]:
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="database insertions.tsv")
    ap.add_argument("--inserts", required=True, help="insertions_inserts.fna")
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-list", default=None,
                    help="scan only these genomes (one basename per line)")
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-ident", type=float, default=0.95)
    ap.add_argument("--min-cov", type=float, default=0.90)
    ap.add_argument("--min-tsd", type=int, default=1)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--batch", type=int, default=100,
                    help="genomes concatenated into one minimap2 reference")
    ap.add_argument("--minimap2", default="minimap2")
    a = ap.parse_args()

    rows = {r["db_id"]: r for r in csv.DictReader(open(a.db), delimiter="\t")
            if r.get("species") == a.species}
    if not rows:
        log("no rows for species %s" % a.species)
        return 0
    ins_all = read_fa(a.inserts)
    ins = {k: v for k, v in ins_all.items() if k in rows}
    log("%s: %d events, %d insert sequences" % (a.species, len(rows), len(ins)))

    # SHARDING. Runtime scales as genomes x inserts: E. faecium (1,038 x
    # 11,087) took ~5h, so E. coli (13,027 x 31,361) is ~175h and does not fit
    # any wall. --genome-list lets one species be split across jobs; each
    # writes its own output and the v5 builder dedups on (bag, flank)
    # globally, so cross-shard duplicates collapse there rather than needing
    # coordination here.
    if a.genome_list:
        genomes = [os.path.basename(x.strip()) for x in open(a.genome_list)
                   if x.strip()]
        log("  %d genomes from %s" % (len(genomes), os.path.basename(a.genome_list)))
    else:
        genomes = sorted(f for f in os.listdir(a.genome_dir) if f.endswith(".fna"))
        log("  %d genomes to scan" % len(genomes))

    # MAPPING DIRECTION. The pilot indexed each genome and mapped all inserts
    # into it: 382 index builds for S. pneumoniae, 24 min. That scales as
    # genomes x inserts and would be ~100h for E. coli. Inverting it -- build
    # ONE .mmi of the inserts, then stream each genome against it -- pays the
    # index cost once. The element's position then comes from the QUERY
    # (contig) coordinates and its coverage from the TARGET (insert) side.
    tmp = tempfile.mkdtemp(prefix="copysite_")
    qf = os.path.join(tmp, "inserts.fna")
    with open(qf, "w") as fh:
        for k, v in ins.items():
            fh.write(">%s\n%s\n" % (k, v))


    st = collections.Counter()
    seen_flank = set()
    out = open(a.out + "_sites.tsv", "w")
    # flank      canonical (min of itself and its revcomp) -- what a bag uses
    # flank_fwd  the SAME window in genomic orientation -- what the round-trip
    #            check needs, since canonicalisation destroys the mapping back
    #            to the contig
    cols = ["db_id", "source_event", "species", "carrier", "contig",
            "elem_start", "elem_end", "tsd_len", "pct_ident", "q_cov",
            "flank", "flank_fwd", "orient", "elem_len",
            "insert_md5", "inserted_len", "empty_site_source"]
    out.write("\t".join(cols) + "\n")

    # BATCHED REFERENCE. Indexing one genome per minimap2 call cost 24 min for
    # 382 genomes, which scales as genomes x inserts and is ~100h for E. coli.
    # Concatenating N genomes into one reference pays the index cost once per
    # batch instead of once per genome, while keeping the inserts as queries.
    # Contig names are prefixed with the carrier so the hit stays attributable.
    batches = [genomes[i:i + a.batch] for i in range(0, len(genomes), a.batch)]
    log("  %d genomes in %d batches of %d" % (len(genomes), len(batches), a.batch))
    for bi, batch in enumerate(batches):
        ref = os.path.join(tmp, "ref.fna")
        ctg = {}
        with open(ref, "w") as fh:
            for g in batch:
                car = g[:-4]
                for name, seq in read_fa(os.path.join(a.genome_dir, g)).items():
                    key = "%s|%s" % (car, name)
                    ctg[key] = seq
                    fh.write(">%s\n%s\n" % (key, seq))
        paf = os.path.join(tmp, "hit.paf")
        # -N/-p OVERRIDE. minimap2 defaults to at most 5 secondary alignments
        # at >=80% of the primary score. With ONE genome per reference that
        # never bound; with 100 genomes an insert hits ~100 homologous loci
        # and the cap silently truncates. Measured: batching without this
        # returned 310,354 hits and 7,441 sites against 1,792,885 and 13,990
        # unbatched -- a 47% loss that looks like a real answer.
        # Finding every copy is the entire point, so secondaries are kept.
        rc = subprocess.run([a.minimap2, "-c", "-x", "asm10", "-t",
                             str(a.threads), "-N", "5000", "-p", "0.1",
                             "--secondary=yes", ref, qf],
                            stdout=open(paf, "w"),
                            stderr=subprocess.DEVNULL).returncode
        if rc != 0:
            st["minimap2_fail"] += 1
            continue
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            # inserts are the QUERIES. Measured: inverting this (genome as
            # query, inserts as reference) returned 7.7x FEWER sites and ran
            # 3.4x slower -- minimap2 reports the best target per query
            # region, so an insert matching a region another insert already
            # claimed is dropped. "Where does each insert occur" needs each
            # insert to be its own query.
            q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            t, ts, te = f[5], int(f[7]), int(f[8])
            nmatch, blen = int(f[9]), int(f[10])
            st["hits"] += 1
            if q not in rows:
                continue
            if blen == 0 or nmatch / blen < a.min_ident:
                st["low_ident"] += 1
                continue
            if (qe - qs) / max(1, qlen) < a.min_cov:
                st["low_cov"] += 1
                continue
            carrier, _, contig = t.partition("|")
            seq = ctg.get(t, "")
            if not seq or te > len(seq):
                st["contig_missing"] += 1
                continue
            if ts < HALF or len(seq) - te < HALF:
                st["too_near_contig_end"] += 1
                continue
            d = tsd_len(seq, ts, te)
            if d < a.min_tsd:
                st["no_tsd"] += 1
                continue
            rec = seq[:ts] + seq[te:]
            if ts < HALF or len(rec) - ts < HALF:
                st["window_would_not_fit"] += 1
                continue
            flank = rec[ts - HALF:ts + HALF]
            if flank.count("N") > 6:
                st["too_many_N"] += 1
                continue
            rcf = LI.rc(flank)
            cf = flank if flank <= rcf else rcf
            orient = "fwd" if flank <= rcf else "rc"
            if cf in seen_flank:
                st["duplicate_flank"] += 1
                continue
            seen_flank.add(cf)
            st["EMITTED"] += 1
            src = rows[q]
            out.write("\t".join(str(x) for x in [
                "%s.copy%d" % (q, st["EMITTED"]), q, a.species, carrier, contig,
                ts, te, d, "%.4f" % (nmatch / blen),
                "%.4f" % ((qe - qs) / max(1, qlen)), cf, flank, orient,
                te - ts,
                src.get("insert_md5", "."), src.get("inserted_len", "."),
                "tsd_reconstructed"]) + "\n")
        log("  batch %d/%d (%d genomes), emitted %d so far"
            % (bi + 1, len(batches), len(batch), st["EMITTED"]))
    out.close()
    log("=" * 70)
    log("%s  emitted %d reconstructed sites" % (a.species, st["EMITTED"]))
    for k in ("hits", "low_ident", "low_cov", "no_tsd", "too_near_contig_end",
              "window_would_not_fit", "too_many_N", "duplicate_flank",
              "contig_missing", "minimap2_fail"):
        if st[k]:
            log("  %-24s %8d" % (k, st[k]))
    log("wrote %s_sites.tsv" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
