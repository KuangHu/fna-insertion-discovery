#!/usr/bin/env python3
"""Can Layer 2 anchors still place at greater outgroup distance? Measure, don't assume.

Layer 3 has never used an outgroup below ~98.2% ANI. That was read as "no distant
genomes available", but the pool already spans the range -- 3,054 genomes sit at
96-97% ANI to cl0000. Two compounding causes kept them out:

  * `80_outgroup_selector.py` ranks candidates by ALIGNED FRACTION descending,
    and aligned fraction falls monotonically with ANI (86% at 99+, 74% at
    96-97%), so the 98-99% band always fills the candidate list first;
  * `--min-af 80` then eliminates almost the whole distant band: of 3,054
    genomes at 96-97%, only 51 pass; at 95-96%, none do.

But whole-genome aligned fraction is the wrong predictor here. It says what
share of the genome aligns end to end; Layer 3 only needs ONE 500 bp anchor pair
at ONE locus to place unambiguously, and anchors sit in conserved flanking
sequence. A genome with 74% AF may still place 90% of anchors.

So this sweep measures the quantity that actually gates Layer 3 -- per ANI band,
the fraction of loci whose anchor pair places uniquely -- using the LOCKED
`place_locus()` so the answer is about the real code path.

Interpretation: if placement holds at 96-97%, the outgroup deficit is a ranking
bug and is free to fix. If placement collapses, the distant band is unusable and
Pool O must be bought at a distance where flanks still align.

    88_anchor_alignability_sweep.py --loci loci.tsv --seed-loci seed.tsv \\
        --ingroup cl0000.manifest --dist ingroup_vs_pool.tsv \\
        --pool-list pool_all.txt --out sweep
"""
import argparse
import csv
import importlib.util
import os
import random
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, write_fasta                        # noqa: E402
from lib.util import log                                        # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_ar", os.path.join(_HERE, "70_allele_reconstructor.py"))
_ar = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ar)

BANDS = [(99.0, 100.1, "99.0+"), (98.0, 99.0, "98.0-99.0"),
         (97.0, 98.0, "97.0-98.0"), (96.0, 97.0, "96.0-97.0"),
         (95.0, 96.0, "95.0-96.0"), (0.0, 95.0, "<95.0")]


def sn(p):
    b = os.path.basename(p)
    for s in (".gz", ".fna", ".fa", ".fasta"):
        if b.endswith(s):
            b = b[: -len(s)]
    return b


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--seed-loci", required=True)
    ap.add_argument("--ingroup", required=True)
    ap.add_argument("--dist", required=True, help="ingroup_vs_pool.tsv")
    ap.add_argument("--pool-list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-band", type=int, default=12)
    ap.add_argument("--max-loci", type=int, default=60)
    ap.add_argument("--anchor", type=int, default=500)
    ap.add_argument("--max-interval", type=int, default=7000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args()

    random.seed(args.seed)
    ing = {sn(l.strip()): l.strip()
           for l in open(args.ingroup) if l.strip()}
    pool = {sn(l.strip()): l.strip()
            for l in open(args.pool_list) if l.strip()}

    def core(n):
        p = n.split("_")
        return p[1].split(".")[0] if len(p) > 1 else n
    icore = {core(x) for x in ing}

    best = {}
    with open(args.dist) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, t = sn(r["Query_file"]), sn(r["Ref_file"])
            if q not in ing or t in ing or core(t) in icore:
                continue
            try:
                a = float(r["ANI"])
                af = min(float(r["Align_fraction_ref"]),
                         float(r["Align_fraction_query"]))
            except (ValueError, KeyError):
                continue
            if t not in best or a > best[t][0]:
                best[t] = (a, af)

    picked = {}
    for lo, hi, name in BANDS:
        cand = [t for t, (a, af) in best.items() if lo <= a < hi and t in pool]
        random.shuffle(cand)
        picked[name] = cand[: args.per_band]
        log("band %-11s available=%-6d sampled=%d" % (name, len(cand), len(picked[name])))

    loci = {r["locus_id"]: r for r in
            csv.DictReader(open(args.loci), delimiter="\t")}
    seeds = {r["locus_id"]: r for r in
             csv.DictReader(open(args.seed_loci), delimiter="\t")}
    lids = [l for l in sorted(loci) if l in seeds][: args.max_loci]
    log("%d loci in the sweep" % len(lids))

    tmp = tempfile.mkdtemp(prefix="sweep.")
    qp = os.path.join(tmp, "anchors.fna")
    fas, kept = {}, []
    with open(qp, "w") as qf:
        for lid in lids:
            sd = seeds[lid]
            gp = ing.get(sd["seed_genome"])
            if not gp:
                continue
            fa = fas.setdefault(sd["seed_genome"], Fasta(gp))
            if sd["contig"] not in fa:
                continue
            a, b = int(sd["start"]), int(sd["end"])
            if a - args.anchor < 0 or b + args.anchor > fa.length(sd["contig"]):
                continue
            write_fasta(qf, lid + "#L", fa.fetch(sd["contig"], a - args.anchor, a))
            write_fasta(qf, lid + "#R", fa.fetch(sd["contig"], b, b + args.anchor))
            kept.append(lid)
    log("%d loci with full anchors" % len(kept))

    rows = []
    for lo, hi, name in BANDS:
        for g in picked[name]:
            paf = os.path.join(tmp, g + ".paf")
            try:
                _ar.map_anchors(qp, pool[g], paf, args.threads, args.minimap2)
            except subprocess.CalledProcessError:
                continue
            hits = _ar.best_hits(paf, 90.0, 0.8)
            uniq = amb = none = 0
            for lid in kept:
                lh, rh = hits.get(lid + "#L", []), hits.get(lid + "#R", [])
                if not lh or not rh:
                    none += 1
                    continue
                pl = _ar.place_locus(lh, rh, args.max_interval, None, 20000, 0.05)
                if pl is None:
                    none += 1
                elif pl["status"] == "ambiguous":
                    amb += 1
                else:
                    uniq += 1
            ani, af = best[g]
            rows.append([name, g, round(ani, 3), round(af, 1), uniq, amb, none,
                         round(100.0 * uniq / max(1, len(kept)), 1)])
            log("  %-11s %-20s ANI=%.2f AF=%.1f  unique=%d amb=%d none=%d (%.0f%%)"
                % (name, g, ani, af, uniq, amb, none,
                   100.0 * uniq / max(1, len(kept))))

    with open(args.out + "_per_genome.tsv", "w") as fh:
        fh.write("band\tgenome\tani\taligned_fraction\tloci_unique\tloci_ambiguous\t"
                 "loci_unplaced\tunique_pct\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    log("=" * 74)
    log("ANCHOR PLACEMENT vs OUTGROUP DISTANCE   (%d loci, locked place_locus)"
        % len(kept))
    log("")
    log("  %-12s %8s %9s %14s %12s" %
        ("ANI band", "genomes", "med AF", "med unique %", "usable?"))
    for lo, hi, name in BANDS:
        v = [r for r in rows if r[0] == name]
        if not v:
            continue
        import statistics
        mu = statistics.median(r[7] for r in v)
        maf = statistics.median(r[3] for r in v)
        verdict = ("yes" if mu >= 70 else "marginal" if mu >= 40 else "no")
        log("  %-12s %8d %9.1f %13.1f%% %12s" % (name, len(v), maf, mu, verdict))
    log("")
    log("  A band is usable for Layer 3 if anchors place uniquely there. Whole-")
    log("  genome aligned fraction is NOT the criterion and should not be the")
    log("  selector's ranking key if these two columns disagree.")
    log("=" * 74)
    log("wrote %s_per_genome.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
