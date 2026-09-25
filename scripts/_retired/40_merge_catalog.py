#!/usr/bin/env python3
"""Stage 40 -- union the two arms into one catalogue. Still no annotation.

    within-genome multi-copy (Arm A)      between-genome empty/filled (Arm B)
                     \\                        /
                      \\                      /
                    universal element catalogue
                                |
                        annotation LAST (stage 50)

The arms have complementary blind spots, which is exactly why both are kept:
  Arm A misses singletons, ancient/diverged copies, and anything the assembler
        collapsed -- and it cannot see the empty allele at all.
  Arm B misses elements absent from every genome in the cluster, and needs a
        cluster-mate that actually lacks the element.
An element found by BOTH is the strongest class: it is repeated within a
genome AND polymorphic between genomes.

POLARITY. Arm B gives presence/absence, not gain/loss. With two assemblies you
cannot tell insertion from deletion -- that is an information limit, not a
software one. So `allele_state` stays PRESENCE_ABSENCE unless --outgroup names
genomes outside the cluster: then the outgroup's allele is taken as ancestral
and the call becomes GAIN or LOSS by parsimony.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log, run                                   # noqa: E402


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(hdr):
                yield dict(zip(hdr, f))


def cluster_sequences(fastas, outdir, mmseqs, min_id, cov, threads):
    """One sequence-identity clustering across BOTH arms, so the same element
    discovered two different ways collapses to one catalogue entry."""
    comb = os.path.join(outdir, "all_elements.fna")
    n = 0
    with open(comb, "w") as out:
        for fa in fastas:
            if not (fa and os.path.exists(fa) and os.path.getsize(fa)):
                continue
            with open(fa) as fh:
                for line in fh:
                    out.write(line)
                    n += line.startswith(">")
    if n == 0:
        return {}, comb
    pref = os.path.join(outdir, "mmseqs")
    tmp = os.path.join(outdir, "mmseqs_tmp")
    run([mmseqs, "easy-linclust", comb, pref, tmp,
         "--min-seq-id", str(min_id), "-c", str(cov),
         "--cov-mode", "1", "--threads", str(threads)])
    rep_of = {}
    ct = pref + "_cluster.tsv"
    if os.path.exists(ct):
        with open(ct) as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")[:2]
                rep_of[mem] = rep
    return rep_of, comb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--armA-glob", default=None,
                    help="glob for *_families.tsv from stage 20")
    ap.add_argument("--armA-elements-glob", default=None,
                    help="glob for *_elements.fna from stage 20")
    ap.add_argument("--armB-refined", default=None, help="armB_refined.tsv")
    ap.add_argument("--armB-events", default=None, help="armB_events.tsv (for carriers)")
    ap.add_argument("--armB-inserts", default=None, help="armB_inserts.fna")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outgroup", default=None,
                    help="file of genome names outside the clusters; enables polarity")
    ap.add_argument("--min-seq-id", type=float, default=0.85)
    ap.add_argument("--cov", type=float, default=0.80)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--mmseqs", default="mmseqs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    outg = set()
    if args.outgroup and os.path.exists(args.outgroup):
        outg = {l.strip() for l in open(args.outgroup) if l.strip()}

    # ---- collect Arm A families
    armA = []
    for path in sorted(glob.glob(args.armA_glob or "")):
        armA.extend(read_tsv(path))
    log("Arm A: %d families" % len(armA))

    # ---- collect Arm B events. Stage 31 is OPTIONAL: stage 30 already emits
    # the inserted sequence and TSD from the graph alleles, so when no refined
    # table exists the events table is promoted to the same schema.
    armB = list(read_tsv(args.armB_refined)) if args.armB_refined \
        and os.path.exists(args.armB_refined) else []
    if not armB and args.armB_events and os.path.exists(args.armB_events):
        log("no stage-31 table: using stage-30 events (graph-resolution coords)")
        for r in read_tsv(args.armB_events):
            fg = (r.get("filled_genomes") or "").split(",")[0]
            eg = (r.get("empty_genomes") or "").split(",")[0]
            armB.append({
                "event_id": r["event_id"], "cluster_id": r["cluster_id"],
                "filled_genome": fg, "empty_genome": eg,
                "empty_chrom": r["chrom"], "empty_pos": r["bubble_start"],
                "insert_len": r.get("insert_len_direct") or r["insert_size"],
                "tsd_len": r.get("tsd_len", "0"),
                "tsd_conf": r.get("tsd_conf", "none"),
                "left_flank": ".", "right_flank": ".",
                "refine_status": "GRAPH_ONLY",
            })
    carriers = {}
    if args.armB_events and os.path.exists(args.armB_events):
        for r in read_tsv(args.armB_events):
            carriers[r["event_id"]] = r
    log("Arm B: %d refined events" % len(armB))

    # ---- one clustering across both arms
    fastas = sorted(glob.glob(args.armA_elements_glob or "")) + \
        ([args.armB_inserts] if args.armB_inserts else [])
    rep_of, _ = cluster_sequences(fastas, args.outdir, args.mmseqs,
                                  args.min_seq_id, args.cov, args.threads)
    log("clustered into %d representatives" % len(set(rep_of.values())))

    # ---- element-level catalogue
    ele = defaultdict(lambda: {"armA": [], "armB": [], "lens": []})
    for r in armA:
        rep = rep_of.get(r["family_id"], r["family_id"])
        ele[rep]["armA"].append(r)
        ele[rep]["lens"].append(int(r["elem_len_consensus"] or r["elem_len_median"]))
    for r in armB:
        if r.get("refine_status") not in ("OK", "GRAPH_ONLY"):
            continue
        rep = rep_of.get(r["event_id"], r["event_id"])
        ele[rep]["armB"].append(r)
        try:
            ele[rep]["lens"].append(int(r["insert_len"]))
        except ValueError:
            pass

    ecols = ["element_id", "rep_seq_id", "len_median", "n_armA_families",
             "n_armB_events", "support", "max_copies_in_one_genome",
             "n_genomes_with_copy", "presence_count", "absence_count",
             "best_S_mobility", "max_tsd_len", "best_tsd_conf",
             "allele_state", "gain_loss_confidence"]
    ef = open(os.path.join(args.outdir, "universal_elements.tsv"), "w")
    ef.write("\t".join(ecols) + "\n")

    vcols = ["event_id", "element_id", "arm", "genome", "cluster_id",
             "chrom", "start", "end", "insert_len", "tsd_len", "tsd_conf",
             "left_flank_200", "right_flank_200", "qc_tier", "allele_state"]
    vf = open(os.path.join(args.outdir, "universal_events.tsv"), "w")
    vf.write("\t".join(vcols) + "\n")

    for i, (rep, d) in enumerate(sorted(ele.items(), key=lambda kv: -len(kv[1]["armA"]) - len(kv[1]["armB"]))):
        eid = "E%05d" % i
        lens = d["lens"] or [0]
        lens.sort()
        med = lens[len(lens) // 2]
        support = ("BOTH_ARMS" if d["armA"] and d["armB"]
                   else "ARM_A_MULTICOPY" if d["armA"] else "ARM_B_EMPTY_FILLED")
        max_copies = max([int(r["n_distinct_loci"]) for r in d["armA"]], default=0)
        genomes = {r["sample_id"] for r in d["armA"]}
        best_S = max([float(r["S_mobility"]) for r in d["armA"]], default=float("nan"))

        pres = absn = 0
        tsd_max, tsd_conf = 0, "none"
        for r in d["armB"]:
            c = carriers.get(r["event_id"])
            if c:
                pres += int(c["n_filled"])
                absn += int(c["n_empty"])
            try:
                if int(r["tsd_len"]) > tsd_max:
                    tsd_max, tsd_conf = int(r["tsd_len"]), r.get("tsd_conf", "none")
            except ValueError:
                pass
            genomes.add(r["filled_genome"])

        # polarity only with an outgroup
        state, conf = "PRESENCE_ABSENCE", "unpolarized"
        if outg and d["armB"]:
            og_filled = sum(1 for r in d["armB"] if r["filled_genome"] in outg)
            og_empty = sum(1 for r in d["armB"] if r["empty_genome"] in outg)
            if og_empty and not og_filled:
                state, conf = "GAIN", "parsimony_outgroup_empty"
            elif og_filled and not og_empty:
                state, conf = "LOSS", "parsimony_outgroup_filled"

        ef.write("\t".join(map(str, [
            eid, rep, med, len(d["armA"]), len(d["armB"]), support,
            max_copies, len(genomes), pres, absn,
            "%.4f" % best_S if best_S == best_S else "NA",
            tsd_max, tsd_conf, state, conf])) + "\n")

        for r in d["armA"]:
            vf.write("\t".join([r["family_id"], eid, "A", r["sample_id"], ".",
                                ".", ".", ".", r["elem_len_consensus"], ".", ".",
                                ".", ".", ".", "MULTICOPY_PRESENT"]) + "\n")
        for r in d["armB"]:
            vf.write("\t".join([r["event_id"], eid, "B", r["filled_genome"],
                                r["cluster_id"], r["empty_chrom"], r["empty_pos"],
                                r["empty_pos"], r["insert_len"], r["tsd_len"],
                                r.get("tsd_conf", "."), r["left_flank"],
                                r["right_flank"], ".", state]) + "\n")
    ef.close()
    vf.close()
    log("wrote universal_elements.tsv (%d elements) + universal_events.tsv" % len(ele))


if __name__ == "__main__":
    main()
