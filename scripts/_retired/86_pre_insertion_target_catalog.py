#!/usr/bin/env python3
"""Layer 4 -- the qualified pre-insertion target catalogue. P1 polarity only.

This module COMBINES; it does not re-derive. Layer 2 already established what
the alleles are and how one becomes the other, and Layer 3 established which way
round it happened. Layer 4 writes those two facts down together:

    Layer 2:  40 bp allele  ->  37 bp + [1279 bp] + 3 bp
    Layer 3:  the 40 bp allele is ancestral
    Layer 4:  pre-insertion target = that exact 40 bp
              inserted sequence    = that exact 1279 bp
              insertion offset     = 37

No new boundary heuristic is introduced here. If a number is not already in the
Layer 2 reconstruction or the Layer 3 polarity call, it does not appear.

TWO CONSISTENCY GATES, both of which refuse rather than repair:

1. POLARITY / DECOMPOSITION COMPATIBILITY. Layer 2 explains the LONG allele as
   the short one plus an insert. If Layer 3 says the LONG allele is ancestral,
   then the event is a loss, and the existing decomposition is not a
   pre-insertion target -- it is the same arithmetic read backwards. Silently
   relabelling it would invent a deletion the pipeline never validated. Such
   loci are emitted as `deletion_candidate` with no target sequence.
   P1 qualifies ancestry; it does not promise the change was an insertion.

2. ALLELE IDENTITY. The ancestral allele arrives from Layer 3 as an ID. IDs are
   assigned per locus by carrier count, and Layer 3b matched outgroups to them
   by sequence similarity, so an ID is a bookkeeping handle rather than a
   sequence. Before any sequence is written out, the ID is resolved back to the
   Layer 2 allele table and its md5 checked against the FASTA. A mismatch is
   fatal for that locus, not a warning.

Coverage note: 72 P1 loci is NOT 72 usable targets. The downstream-ready
denominator is

    P1 polarity  AND  resolved decomposition  AND  usable target context

and it is reported separately, because that is the number an RNA-guide dataset
can actually be built from.

    86_pre_insertion_target_catalog.py --polarity polarity/all_cl0000 \\
        --recon allele_recon/armB_cl0000 --out catalog/cl0000
"""
import argparse
import csv
import hashlib
import os
import statistics
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, revcomp                            # noqa: E402
from lib.util import log                                        # noqa: E402


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()[:12] if s else "EMPTY"


def read_fasta_alleles(path):
    out, cur = defaultdict(dict), None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur:
                out[cur[0]][cur[1]] = []
        elif cur:
            out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "") for a, v in d.items()}
            for k, d in out.items()}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polarity", nargs="+", required=True,
                    help="Layer 3 output dirs (polarity.tsv)")
    ap.add_argument("--recon", nargs="+", required=True,
                    help="Layer 2 dirs, same order as --polarity")
    ap.add_argument("--refinements", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--genomes", nargs="*", default=[],
                    help="manifests of the carrier genomes. With these, context "
                         "is RE-EXTRACTED from the source FASTA centred on the "
                         "junction, instead of being cut out of the padded "
                         "anchored interval. Recovers loci whose window would "
                         "otherwise truncate at the anchor, without relaxing "
                         "any quality criterion.")
    ap.add_argument("--context-windows", default="50,100,200,500")
    ap.add_argument("--tiers", default="P1",
                    help="polarity tiers to emit targets for (default P1 only)")
    args = ap.parse_args()

    keep = set(args.tiers.split(","))
    os.makedirs(args.out, exist_ok=True)

    pol, loci, atab, aseq, ctx = {}, {}, defaultdict(dict), {}, defaultdict(dict)
    cluster_of = {}
    for pdir, rdir in zip(args.polarity, args.recon):
        cid = os.path.basename(rdir).replace("armB_", "")
        p = os.path.join(pdir, "polarity.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                pol[r["locus_id"]] = r
                cluster_of[r["locus_id"]] = cid
        p = os.path.join(rdir, "loci.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                loci[r["locus_id"]] = r
        p = os.path.join(rdir, "alleles.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                atab[r["locus_id"]][r["allele_id"]] = r
        aseq.update(read_fasta_alleles(os.path.join(rdir, "alleles.fna")))
        p = os.path.join(rdir, "target_context.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                ctx[(r["locus_id"], r["allele_id"])][int(r["window_bp"])] = r

    windows = [int(x) for x in args.context_windows.split(",") if x.strip()]
    gpaths = {}
    for m in args.genomes:
        if not os.path.exists(m):
            continue
        for line in open(m):
            line = line.strip()
            if line:
                b = os.path.basename(line)
                for suf in (".gz", ".fna", ".fa", ".fasta"):
                    if b.endswith(suf):
                        b = b[: -len(suf)]
                gpaths[b] = line
    fas = {}

    def genome(g):
        if g not in fas:
            fas[g] = Fasta(gpaths[g]) if g in gpaths else None
        return fas[g]

    focal = {}
    for f in args.refinements:
        if not os.path.exists(f):
            continue
        cid = os.path.basename(os.path.dirname(f))
        for r in csv.DictReader(open(f), delimiter="\t"):
            focal[cid] = r.get("focal_members", "")

    cand = [l for l, r in pol.items() if r.get("polarity_tier") in keep]
    log("%d loci at tier(s) %s" % (len(cand), ",".join(sorted(keep))))

    rows, rej, cls = [], Counter(), Counter()
    for lid in sorted(cand):
        P, L = pol[lid], loci.get(lid)
        if L is None:
            rej["no_layer2_record"] += 1
            continue
        anc = P.get("ancestral_allele", ".")
        short_id, long_id = L.get("shortest_allele"), L.get("longest_allele")

        # ---- gate 2: the ID must resolve to a real, verified sequence -------
        arec = atab.get(lid, {}).get(anc)
        aseq_ = aseq.get(lid, {}).get(anc)
        if arec is None or aseq_ is None:
            rej["ancestral_allele_id_unresolvable"] += 1
            continue
        if arec.get("md5") != md5(aseq_) or int(arec["length"]) != len(aseq_):
            rej["allele_hash_mismatch"] += 1
            continue

        # ---- gate 1: does the polarity match the decomposition direction? ---
        if anc == long_id and short_id != long_id:
            # ancestral is the LONG allele -> the event removed sequence.
            # Layer 2 never validated a deletion model, so no target is emitted.
            cls["deletion_candidate"] += 1
            rows.append({
                "event_id": lid, "discovery_cluster": cluster_of.get(lid, "."),
                "focal_clade": focal.get(cluster_of.get(lid, ""), "."),
                "pre_event_allele_id": anc, "derived_allele_id": short_id,
                "polarity_tier": P.get("polarity_tier"),
                "polarity_support": P.get("ancestral_bootstrap_support"),
                "sampling_stability": "not_tested",
                "event_direction": "deletion_candidate",
                "pre_event_interval_len": len(aseq_),
                "downstream_ready": "no",
                "reject_reason": "ancestral_is_long_allele_no_insertion_model"})
            continue
        if anc != short_id:
            rej["ancestral_allele_not_in_decomposition"] += 1
            cls["complex_polarity"] += 1
            continue

        ev = L.get("event_class", "")
        if ev not in ("insertion_target_retained", "replacement", "pure_insertion"):
            rej["decomposition_not_an_insertion_model"] += 1
            continue

        derived_seq = aseq.get(lid, {}).get(long_id, "")
        lcp, lcs = int(L.get("lcp_bp", 0) or 0), int(L.get("lcs_bp", 0) or 0)
        ins_len = int(L.get("inserted_len", 0) or 0)
        ins_seq = derived_seq[lcp:len(derived_seq) - lcs] if derived_seq else ""
        if ins_seq and len(ins_seq) != ins_len:
            rej["inserted_length_inconsistent"] += 1
            continue

        # Context must be centred on the INSERTION POINT, not on the anchored
        # interval. The loci were seeded with a 200 bp pad each side, so the
        # interval between anchors is ~409 bp of ordinary flanking sequence
        # with the junction sitting at offset lcp inside it. Layer 2's
        # target_context.tsv brackets the whole interval, which for an
        # RNA-guide comparison would centre the window ~200 bp away from the
        # actual target. lcp is already a Layer 2 quantity, so this re-centres
        # rather than introducing a new boundary estimate.
        pt = lcp
        winc, trunc, src = {}, False, "anchored_interval"
        # Preferred route: go back to the carrier genome. Layer 2 already
        # recorded which genome carries the pre-event allele, on which contig
        # and strand, and the junction's absolute coordinate, so this reads more
        # sequence around a known point -- it does not re-infer the point.
        cg = L.get("pre_event_genome", ".")
        cc = L.get("pre_event_contig", ".")
        cs = L.get("pre_event_strand", "+")
        try:
            jabs = int(L.get("junction_L_abs", -1))
        except ValueError:
            jabs = -1
        fa = genome(cg) if cg not in (".", "") else None
        if fa is not None and cc in fa and jabs >= 0:
            clen = fa.length(cc)
            for w in windows:
                lo, hi = jabs - w, jabs + w
                if lo < 0 or hi > clen:
                    trunc = True
                    lo, hi = max(0, lo), min(clen, hi)
                lseq, rseq = fa.fetch(cc, lo, jabs), fa.fetch(cc, jabs, hi)
                if cs == "-":
                    lseq, rseq = revcomp(rseq), revcomp(lseq)
                winc[w] = (lseq, rseq)
            src = "carrier_genome"
        else:
            for w in windows:
                lo, hi = pt - w, pt + w
                if lo < 0 or hi > len(aseq_):
                    trunc = True
                winc[w] = (aseq_[max(0, lo):pt], aseq_[pt:min(len(aseq_), hi)])
        # the target interval proper: what sits between the two junctions.
        # For a clean insertion with retention this is empty and the event is a
        # point insertion with `ambiguity` bp of microhomology.
        tgt = aseq_[lcp:len(aseq_) - lcs] if len(aseq_) - lcs > lcp else ""
        usable = not trunc
        # END-TO-END CHECK. If the ancestry and the insertion are both right,
        # splicing the insert into the pre-event allele at the inferred offset
        # must regenerate the observed derived allele. This is the one test that
        # exercises Layer 2 and Layer 3 together on the final product.
        rebuilt = aseq_[:lcp] + ins_seq + aseq_[len(aseq_) - lcs:] \
            if lcs else aseq_[:lcp] + ins_seq
        rec_ident = rec_cov = 0.0
        if derived_seq:
            n_match = sum(1 for a, b in zip(rebuilt, derived_seq) if a == b)
            rec_ident = n_match / max(1, len(derived_seq))
            rec_cov = min(len(rebuilt), len(derived_seq)) / max(
                1, max(len(rebuilt), len(derived_seq)))
        cls[ev] += 1
        rows.append({
            "event_id": lid, "discovery_cluster": cluster_of.get(lid, "."),
            "focal_clade": focal.get(cluster_of.get(lid, ""), "."),
            "pre_event_allele_id": anc, "derived_allele_id": long_id,
            "polarity_tier": P.get("polarity_tier"),
            "polarity_support": P.get("ancestral_bootstrap_support"),
            "sampling_stability": "no_flip_observed",
            "event_direction": "insertion",
            "event_class": ev,
            "pre_event_interval_seq": aseq_,
            "pre_event_interval_len": len(aseq_),
            "inserted_seq": ins_seq, "inserted_len": len(ins_seq),
            "insertion_offset": lcp,
            "target_bases_retained_left": lcp,
            "target_bases_retained_right": lcs,
            "target_bases_lost": L.get("target_bases_lost", "0"),
            "insertion_point_offset": pt,
            "target_interval_seq": tgt or "-",
            "target_interval_len": len(tgt),
            "pre_insertion_context_50": winc.get(50, ("", ""))[0] + "|" +
                                        winc.get(50, ("", ""))[1],
            "pre_insertion_context_100": winc.get(100, ("", ""))[0] + "|" +
                                         winc.get(100, ("", ""))[1],
            "pre_insertion_context_200": winc.get(200, ("", ""))[0] + "|" +
                                         winc.get(200, ("", ""))[1],
            "pre_insertion_context_500": winc.get(500, ("", ""))[0] + "|" +
                                         winc.get(500, ("", ""))[1],
            "context_truncated": "yes" if trunc else "no",
            "context_source": src,
            "allele_reconstruction_identity": "%.4f" % rec_ident,
            "allele_reconstruction_coverage": "%.4f" % rec_cov,
            "context_source_genome":
                ctx.get((lid, anc), {}).get(50, {}).get("source_genome", "."),
            "placement_status": L.get("placement_status", "."),
            "junction_ambiguity_bp": L.get("junction_ambiguity_bp", "."),
            "decomposition_quality": L.get("decomposition_quality", "."),
            "decomposition_method": L.get("decomposition_method", "."),
            "validation_flag": ".",
            "downstream_ready": "yes" if (usable and rec_ident >= 0.99) else "no",
            "reject_reason": ("" if (usable and rec_ident >= 0.99)
                              else "context_truncated" if not usable
                              else "allele_reconstruction_failed")})

    cols = ["event_id", "discovery_cluster", "focal_clade", "pre_event_allele_id",
            "derived_allele_id", "polarity_tier", "polarity_support",
            "sampling_stability", "event_direction", "event_class",
            "pre_event_interval_seq", "pre_event_interval_len",
            "insertion_point_offset", "target_interval_seq", "target_interval_len",
            "context_truncated", "inserted_seq",
            "inserted_len", "insertion_offset", "target_bases_retained_left",
            "target_bases_retained_right", "target_bases_lost",
            "pre_insertion_context_50", "pre_insertion_context_100",
            "pre_insertion_context_200", "pre_insertion_context_500",
            "context_source", "allele_reconstruction_identity",
            "allele_reconstruction_coverage", "context_source_genome",
            "placement_status", "junction_ambiguity_bp", "decomposition_quality",
            "decomposition_method", "validation_flag", "downstream_ready",
            "reject_reason"]
    with open(os.path.join(args.out, "pre_insertion_targets.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with open(os.path.join(args.out, "pre_insertion_targets.fna"), "w") as fh:
        for r in rows:
            if r.get("downstream_ready") == "yes":
                fh.write(">%s|pre_event|len=%d\n%s\n"
                         % (r["event_id"], r["pre_event_interval_len"],
                            r["pre_event_interval_seq"]))
    with open(os.path.join(args.out, "inserted_elements.fna"), "w") as fh:
        for r in rows:
            if r.get("downstream_ready") == "yes" and r.get("inserted_seq"):
                fh.write(">%s|inserted|len=%d\n%s\n"
                         % (r["event_id"], r["inserted_len"], r["inserted_seq"]))

    ready = [r for r in rows if r.get("downstream_ready") == "yes"]
    log("=" * 74)
    log("PRE-INSERTION TARGET CATALOGUE")
    log("  P1 loci considered              %4d" % len(cand))
    for k, v in cls.most_common():
        log("    %-30s %4d" % (k, v))
    if rej:
        log("")
        log("  refused:")
        for k, v in rej.most_common():
            log("    %-38s %4d" % (k, v))
    log("")
    log("  DOWNSTREAM-READY  (P1 AND resolved decomposition AND target context)")
    log("    %4d of %d P1 loci  =  %.1f%%"
        % (len(ready), len(cand), 100.0 * len(ready) / len(cand) if cand else 0))
    if ready:
        il = [r["inserted_len"] for r in ready]
        pl = [r["pre_event_interval_len"] for r in ready]
        lost = sum(1 for r in ready if str(r["target_bases_lost"]) == "0")
        amb = sum(1 for r in ready if str(r["junction_ambiguity_bp"]) == "0")
        log("")
        log("    median inserted length        %6d bp  (range %d-%d)"
            % (statistics.median(il), min(il), max(il)))
        log("    median pre-event interval     %6d bp  (range %d-%d)"
            % (statistics.median(pl), min(pl), max(pl)))
        log("    target_bases_lost == 0        %6.1f%%" % (100.0 * lost / len(ready)))
        log("    junction_ambiguity == 0       %6.1f%%" % (100.0 * amb / len(ready)))
        tl = [r["target_interval_len"] for r in ready]
        pt_ = [r["insertion_point_offset"] for r in ready]
        log("    target interval (between junctions): median %d bp, max %d"
            % (statistics.median(tl), max(tl)))
        log("    insertion point offset within interval: median %d bp"
            % statistics.median(pt_))
        ri = [float(r["allele_reconstruction_identity"]) for r in ready]
        log("    allele reconstruction identity: median %.4f, min %.4f"
            % (statistics.median(ri), min(ri)))
        cs_ = Counter(r["context_source"] for r in ready)
        log("    context source: " + ", ".join("%s=%d" % kv for kv in cs_.items()))
        hist = Counter(il)
        log("    inserted-length modes: " + " ".join(
            "%dbp:%d" % (k, v) for k, v in hist.most_common(6)))
    log("=" * 74)
    log("wrote %s/pre_insertion_targets.{tsv,fna} + inserted_elements.fna"
        % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
