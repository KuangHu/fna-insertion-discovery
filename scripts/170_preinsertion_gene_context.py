#!/usr/bin/env python3
"""What was at the target site BEFORE the element arrived -- gene or not?

Every decomposed event records where its EMPTY allele sits in a real genome:
pre_event_genome, pre_event_contig, and junction_L_abs == junction_R_abs, a
single base. So the pre-insertion sequence can be read from an uninterrupted
genome and its genes called there.

WHY THE EMPTY ALLELE AND NOT THE CARRIER. Calling genes on a genome that HAS
the element shows a gene already broken by it -- the answer would be built
into the input. The empty allele is the same locus in a genome that never
received the insertion, so its gene structure is the pre-insertion state. This
is the one question this pipeline is unusually well placed to answer, because
it reconstructs the empty allele rather than inferring it.

WHY A GENOMIC WINDOW AND NOT target_seq. target_seq has median length 409 bp
(it is the anchored interval, mostly the 200 bp pad), far shorter than a gene.
prodigal on it would call fragments and every insertion would look intergenic
or mid-gene by artefact. A window of +/-W kb from the real contig contains
whole genes with real start and stop codons.

CLASSIFICATION of where the junction falls, on the empty sequence:

  intragenic      strictly inside a CDS -- the insertion disrupts a gene
                    [gene]  ->  [gene 5' part][INSERT][gene 3' part]
  intergenic      between two CDS
  gene_boundary   within --edge bp of a CDS start or stop, where prodigal's
                  boundary uncertainty makes intragenic/intergenic a coin flip
  no_gene_called  prodigal found nothing in the window

For intragenic hits the module reports which gene, how far into it the element
landed (as a fraction of the CDS), the strand, and the sizes of the two
fragments the insertion creates -- i.e. the architecture.

WHAT IS NOT CLAIMED. prodigal is ab initio; a called CDS is an open reading
frame, not an annotated or expressed gene. No functional annotation is used,
consistent with the rest of the pipeline.

AND THE INTERPRETATION THAT MUST NOT BE MADE. A low intragenic fraction is NOT
evidence that the element avoids genes. Every event here comes from a cell
that survived: an insertion into an essential gene kills its host and is never
sampled. What is observed is (where elements land) x (whether the cell lived),
and this module cannot separate the two factors. The matched null below shares
the same limit -- it says where an insertion COULD land in that window, so the
gap between null and observed is dominated by selection, not by targeting
preference.

Separating them needs insertions that selection has not yet filtered. The
usable handle is FREQUENCY: a singleton present in one genome is recent, one
fixed across the species is old and has survived. If the intragenic fraction
falls as frequency rises, selection is acting in front of us; if it is flat,
selection is not the explanation. `carriers_derived` and the M3 presence
counts supply that stratification.
"""
import argparse
import collections
import csv
import glob
import os
import random
import subprocess
import sys
import tempfile

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--events", required=True, help="dedup loci_to_event.tsv")
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=5000)
    ap.add_argument("--edge", type=int, default=15,
                    help="bp from a CDS end where the call is ambiguous")
    ap.add_argument("--prodigal", default="prodigal")
    a = ap.parse_args()

    ev_of = {}
    for r in csv.DictReader(open(a.events), delimiter="\t"):
        ev_of[r["locus_id"]] = r["event_id"]

    dirs = []
    for pat in a.recon:
        dirs.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[")
                    else [pat])

    # one row per EVENT, not per locus: the same event seen by several panels
    # is one insertion and must not be counted repeatedly.
    want, seen_ev = [], set()
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import lib_insert as LI
    for d in dirs:
        p = os.path.join(d, "loci.tsv")
        if not os.path.exists(p):
            continue
        panel = LI.panel_label(d)
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["event_class"] not in ("insertion_target_retained", "replacement"):
                continue
            g, c = r.get("pre_event_genome", "."), r.get("pre_event_contig", ".")
            if g in (".", "") or c in (".", ""):
                continue
            try:
                jl, jr = int(r["junction_L_abs"]), int(r["junction_R_abs"])
            except (ValueError, KeyError):
                continue
            lid = LI.qualify(panel, r["locus_id"])
            e = ev_of.get(lid)
            if not e or e in seen_ev:
                continue
            seen_ev.add(e)
            want.append({"event": e, "genome": g, "contig": c,
                         "jl": jl, "jr": jr,
                         "ilen": r.get("inserted_len", "."),
                         "overlap": r.get("junction_ambiguity_bp", "0")})
    log("%s: %d distinct events with a pre-insertion locus" % (a.species, len(want)))
    if not want:
        return 0

    by_genome = collections.defaultdict(list)
    for w in want:
        by_genome[w["genome"]].append(w)
    log("  across %d genomes" % len(by_genome))

    rng = random.Random(17)   # fixed seed: the null must be reproducible
    tmp = tempfile.mkdtemp(prefix="precontext_")
    st = collections.Counter()
    out = open(a.out + "_context.tsv", "w")
    cols = ["event_id", "species", "genome", "contig", "junction_pos",
            "inserted_len", "junction_overlap", "call", "cds_id",
            "cds_start", "cds_end", "cds_len", "cds_strand",
            "pos_in_cds", "frac_into_cds", "frag5_len", "frag3_len",
            "dist_to_nearest_cds"]
    out.write("\t".join(cols) + "\n")

    for gi, (gname, ws) in enumerate(sorted(by_genome.items())):
        gp = os.path.join(a.genome_dir, gname + ".fna")
        if not os.path.exists(gp):
            st["genome_missing"] += len(ws)
            continue
        ctg = read_fa(gp)
        # one prodigal call per genome over all its windows, not per event
        wf = os.path.join(tmp, "w.fna")
        meta = []
        with open(wf, "w") as fh:
            for i, w in enumerate(ws):
                seq = ctg.get(w["contig"], "")
                if not seq:
                    st["contig_missing"] += 1
                    continue
                lo = max(0, w["jl"] - a.window)
                hi = min(len(seq), w["jl"] + a.window)
                sub = seq[lo:hi]
                if len(sub) < 500:
                    st["window_too_small"] += 1
                    continue
                nm = "w%d" % i
                meta.append((nm, w, lo, len(sub)))
                fh.write(">%s\n%s\n" % (nm, sub))
        if not meta:
            continue
        gff = os.path.join(tmp, "w.gff")
        rc = subprocess.run([a.prodigal, "-i", wf, "-f", "gff", "-o", gff,
                             "-p", "meta", "-q"],
                            stderr=subprocess.DEVNULL).returncode
        if rc != 0:
            st["prodigal_fail"] += len(meta)
            continue
        cds = collections.defaultdict(list)
        for line in open(gff):
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[2] != "CDS":
                continue
            cds[f[0]].append((int(f[3]) - 1, int(f[4]), f[6]))
        for nm, w, lo, slen in meta:
            # MATCHED NULL. "33% intragenic" is uninterpretable on its own:
            # bacterial genomes are ~85-90% coding, so random insertion would
            # be mostly intragenic. The null is a position drawn uniformly
            # from the SAME window and classified by the SAME rule, which
            # controls for that window's gene density, prodigal's behaviour on
            # a 10kb fragment, and the edge band all at once. The comparison
            # that means something is observed vs this, not vs a textbook
            # coding fraction.
            hits_n = cds.get(nm, [])
            if hits_n:
                jn = rng.randrange(a.edge, max(a.edge + 1, slen - a.edge))
                ins_n = [(s2, e2) for s2, e2, _ in hits_n if s2 <= jn < e2]
                if ins_n:
                    s2, e2 = max(ins_n, key=lambda x: x[1] - x[0])
                    o2, L2 = jn - s2, e2 - s2
                    st["null_" + ("gene_boundary" if (o2 < a.edge or L2 - o2 < a.edge)
                                  else "intragenic")] += 1
                else:
                    d2n = min(min(abs(jn - s2), abs(jn - e2))
                              for s2, e2, _ in hits_n)
                    st["null_" + ("gene_boundary" if d2n < a.edge
                                  else "intergenic")] += 1
                st["null_checked"] += 1

            # junction position inside the extracted window
            j = w["jl"] - lo
            hits = cds.get(nm, [])
            st["checked"] += 1
            if not hits:
                st["no_gene_called"] += 1
                out.write("\t".join(str(x) for x in [
                    w["event"], a.species, w["genome"], w["contig"], w["jl"],
                    w["ilen"], w["overlap"], "no_gene_called",
                    ".", ".", ".", ".", ".", ".", ".", ".", ".", "."]) + "\n")
                continue
            inside = [(s, e, sd) for s, e, sd in hits if s <= j < e]
            if inside:
                s, e, sd = max(inside, key=lambda x: x[1] - x[0])
                off = j - s
                L = e - s
                near_edge = off < a.edge or (L - off) < a.edge
                call = "gene_boundary" if near_edge else "intragenic"
                st[call] += 1
                # the architecture: the insertion splits the CDS in two
                f5, f3 = off, L - off
                if sd == "-":
                    f5, f3 = f3, f5
                out.write("\t".join(str(x) for x in [
                    w["event"], a.species, w["genome"], w["contig"], w["jl"],
                    w["ilen"], w["overlap"], call,
                    "%s:%d-%d" % (w["contig"], lo + s, lo + e),
                    lo + s, lo + e, L, sd, off, "%.4f" % (off / L),
                    f5, f3, 0]) + "\n")
            else:
                d2 = min(min(abs(j - s), abs(j - e)) for s, e, _ in hits)
                call = "gene_boundary" if d2 < a.edge else "intergenic"
                st[call] += 1
                out.write("\t".join(str(x) for x in [
                    w["event"], a.species, w["genome"], w["contig"], w["jl"],
                    w["ilen"], w["overlap"], call,
                    ".", ".", ".", ".", ".", ".", ".", ".", ".", d2]) + "\n")
        if (gi + 1) % 500 == 0:
            log("  %d/%d genomes" % (gi + 1, len(by_genome)))
    out.close()

    n = max(1, st["checked"])
    log("=" * 70)
    log("PRE-INSERTION GENE CONTEXT -- %s" % a.species)
    log("  events checked        %8d" % st["checked"])
    for k in ("intragenic", "intergenic", "gene_boundary", "no_gene_called"):
        log("    %-18s  %8d  %5.1f%%" % (k, st[k], 100.0 * st[k] / n))
    nn = max(1, st["null_checked"])
    log("")
    log("  MATCHED NULL -- a uniform random position in the SAME windows")
    for k in ("intragenic", "intergenic", "gene_boundary"):
        log("    %-18s  %8d  %5.1f%%   (observed %5.1f%%)"
            % (k, st["null_" + k], 100.0 * st["null_" + k] / nn,
               100.0 * st[k] / n))
    log("")
    for k in ("genome_missing", "contig_missing", "window_too_small",
              "prodigal_fail"):
        if st[k]:
            log("  %-22s %8d" % (k, st[k]))
    log("")
    log("  intragenic = the element landed INSIDE a coding region:")
    log("    [gene]  ->  [gene 5'][INSERT][gene 3']")
    log("  prodigal is ab initio, so a CDS is an ORF, not an annotated gene.")
    log("")
    log("  A low intragenic fraction is NOT evidence the element avoids genes.")
    log("  Every event here survived: insertions into essential genes kill the")
    log("  host and are never sampled. Observed = targeting x survival, and")
    log("  neither this table nor the null separates the two.")
    log("=" * 70)
    log("wrote %s_context.tsv" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
