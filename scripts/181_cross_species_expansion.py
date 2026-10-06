#!/usr/bin/env python3
"""Find every background a candidate element occupies -- across ALL species.

The evidence a candidate system needs is the SAME element in MANY INDEPENDENT
TARGET BACKGROUNDS. Those backgrounds do not have to come from one species, and
waiting for a species to accumulate enough complete genomes is the slow way to
get them. This searches the existing candidates against every assembly we hold,
regardless of which species discovered them.

DISCOVERY STAYS FAMILY-AGNOSTIC. Homology is used here only to EXPAND evidence
for a candidate that structural comparison already found. No element model ever
selects a locus, which is the project's founding constraint.

QUERY SET: one representative insert per CDS cluster (8,672) rather than all
105,048 inserts. Members of a cluster are the same element family, so querying
all of them finds the same loci repeatedly at 12x the cost.

LOCUS-LEVEL QUALIFICATION, not genome-level. A globally fragmented assembly can
still carry a complete left-flank | element | right-flank, and a high N50 does
not guarantee THIS locus is intact. So each hit is judged on its own:

  flank_complete   >= --margin bp of contig on BOTH sides of the element
  truncated_left / truncated_right   the element runs into a contig end, so its
                   extent is unknown and the flank is not evidence
  A MISSING OR BROKEN REGION IS NEVER EMPTY-SITE EVIDENCE. An element at a
  contig edge is recorded as truncated and excluded from empty-site work rather
  than being treated as an absence.

INDEPENDENCE. Backgrounds are counted per (species, genome, locus), and the
same assembly re-deposited under two accessions, or two assemblies of one
BioSample, are not two backgrounds. Genome identity is reported so the caller
can collapse on BioSample; this script does not assume accession == sample.

CALIBRATION DOES NOT TRANSFER. Thresholds tuned on within-species comparison do
not apply to a cross-species evidence set -- identity between two species'
copies of an element is a different quantity from identity within one. The
per-hit identity and coverage are written out so a cross-species score can be
calibrated separately rather than inheriting the within-species cut.
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


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inserts", required=True)
    ap.add_argument("--cds", required=True, help="insert_cds.tsv, db_id -> cluster")
    ap.add_argument("--db", required=True, help="insertions.tsv, for source species")
    ap.add_argument("--target-dir", required=True, help="genomes of ONE target species")
    ap.add_argument("--target-species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--margin", type=int, default=500,
                    help="bp of contig required on BOTH sides to call a flank complete")
    ap.add_argument("--flank", type=int, default=200, help="bp of flank to record")
    ap.add_argument("--min-ident", type=float, default=0.80,
                    help="lower than the within-species 0.95: a cross-species "
                         "copy of the same element is genuinely more diverged")
    ap.add_argument("--min-cov", type=float, default=0.80)
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--minimap2", default="minimap2")
    a = ap.parse_args()

    clu = {}
    for r in csv.DictReader(open(a.cds), delimiter="\t"):
        c = r.get("cds_cluster_id", ".")
        if c and c != ".":
            clu[r["db_id"]] = c
    src_sp = {r["db_id"]: r.get("species", ".")
              for r in csv.DictReader(open(a.db), delimiter="\t")}
    ins = read_fa(a.inserts)

    # one representative per cluster: the longest member, so a partial copy
    # does not become the query that defines the family
    rep = {}
    for dbid, c in clu.items():
        s = ins.get(dbid)
        if not s:
            continue
        if c not in rep or len(s) > len(ins[rep[c]]):
            rep[c] = dbid
    log("%d CDS clusters, %d with a representative insert" % (
        len(set(clu.values())), len(rep)))

    tmp = tempfile.mkdtemp(prefix="xspec_")
    qf = os.path.join(tmp, "reps.fna")
    with open(qf, "w") as fh:
        for c, dbid in sorted(rep.items()):
            fh.write(">%s|%s|%s\n%s\n" % (c, dbid, src_sp.get(dbid, "."), ins[dbid]))

    genomes = sorted(f for f in os.listdir(a.target_dir) if f.endswith(".fna"))
    log("target %s: %d genomes" % (a.target_species, len(genomes)))
    batches = [genomes[i:i + a.batch] for i in range(0, len(genomes), a.batch)]

    st = collections.Counter()
    seen = set()
    out = open(a.out + "_hits.tsv", "w")
    cols = ["cds_cluster", "source_event", "source_species", "target_species",
            "genome", "contig", "elem_start", "elem_end", "elem_len",
            "pct_ident", "q_cov", "cross_species", "flank_complete",
            "dist_left", "dist_right", "left_flank", "right_flank", "tsd_len"]
    out.write("\t".join(cols) + "\n")

    for bi, batch in enumerate(batches):
        ref = os.path.join(tmp, "ref.fna")
        ctg = {}
        with open(ref, "w") as fh:
            for g in batch:
                car = g[:-4]
                for nm, sq in read_fa(os.path.join(a.target_dir, g)).items():
                    k = "%s|%s" % (car, nm)
                    ctg[k] = sq
                    fh.write(">%s\n%s\n" % (k, sq))
        paf = os.path.join(tmp, "h.paf")
        if subprocess.run([a.minimap2, "-c", "-x", "asm20", "-t",
                           str(a.threads), "-N", "5000", "-p", "0.1",
                           "--secondary=yes", ref, qf],
                          stdout=open(paf, "w"),
                          stderr=subprocess.DEVNULL).returncode != 0:
            st["minimap2_fail"] += 1
            continue
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            t, ts, te = f[5], int(f[7]), int(f[8])
            nmatch, blen = int(f[9]), int(f[10])
            st["hits"] += 1
            if blen == 0 or nmatch / blen < a.min_ident:
                st["low_ident"] += 1
                continue
            if (qe - qs) / max(1, qlen) < a.min_cov:
                st["low_cov"] += 1
                continue
            cls, dbid, ssp = q.split("|", 2)
            carrier, _, contig = t.partition("|")
            seq = ctg.get(t, "")
            if not seq or te > len(seq):
                st["contig_missing"] += 1
                continue
            dl, dr = ts, len(seq) - te
            # LOCUS-LEVEL qualification. An element running into a contig end
            # has unknown extent; its flank is not evidence and the absence of
            # sequence beyond it is not an empty site.
            complete = dl >= a.margin and dr >= a.margin
            if not complete:
                st["truncated_left" if dl < a.margin else "truncated_right"] += 1
            lf = seq[max(0, ts - a.flank):ts]
            rf = seq[te:te + a.flank]
            n = 0
            while n < 40 and ts - 1 - n >= 0 and te - 1 - n >= 0 \
                    and seq[ts - 1 - n] == seq[te - 1 - n]:
                n += 1
            key = (cls, carrier, contig, ts // 50)
            if key in seen:
                st["dup_locus"] += 1
                continue
            seen.add(key)
            xs = "yes" if ssp != a.target_species else "no"
            st["EMITTED"] += 1
            st["cross" if xs == "yes" else "same"] += 1
            if complete:
                st["complete_flank"] += 1
            out.write("\t".join(str(x) for x in [
                cls, dbid, ssp, a.target_species, carrier, contig, ts, te,
                te - ts, "%.4f" % (nmatch / blen),
                "%.4f" % ((qe - qs) / max(1, qlen)), xs,
                "yes" if complete else "no", dl, dr, lf, rf, n]) + "\n")
        log("  batch %d/%d  emitted %d" % (bi + 1, len(batches), st["EMITTED"]))
    out.close()

    log("=" * 72)
    log("CROSS-SPECIES EXPANSION -- target %s" % a.target_species)
    log("  hits %d -> emitted %d" % (st["hits"], st["EMITTED"]))
    log("    from a DIFFERENT source species   %8d" % st["cross"])
    log("    from the same species             %8d" % st["same"])
    log("    with a complete flank (>=%dbp)    %8d" % (a.margin, st["complete_flank"]))
    for k in ("low_ident", "low_cov", "truncated_left", "truncated_right",
              "dup_locus", "contig_missing", "minimap2_fail"):
        if st[k]:
            log("  %-22s %8d" % (k, st[k]))
    log("")
    log("  A truncated locus is NOT an empty site. Its element extent is")
    log("  unknown and the missing sequence is absent data, not absence.")
    log("=" * 72)
    log("wrote %s_hits.tsv" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
