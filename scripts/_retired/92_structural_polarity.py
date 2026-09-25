#!/usr/bin/env python3
"""Polarity from target-site duplication alone. Two genomes, no panel, no tree.

Layer 3 answers one question -- which of two alleles is ancestral -- and pays for
it with a per-locus ML genealogy, a monophyletic focal clade and an informative
outgroup panel. That cost is why 234 Layer-2 loci become 20 usable targets.
And measured: **79.6% of P1 loci are constant within the focal clade**, so in
four cases out of five the tree contributes nothing and the answer comes from
outgroup state alone.

For any element that duplicates its target site, the answer is already inside the
allele pair:

    ancestral :  ... P [tsd] S ...
    derived   :  ... P [tsd] ELEMENT [tsd] S ...

One copy becoming two is not reversible by the same mechanism, so the direct
repeat flanking the insert IS the polarity. No outgroup, no clade, no tree.

THIS IS STRUCTURE, NOT ANNOTATION. No HMM, no IS library, no family model; the
repeat is measured inside the two alleles being compared and the test is
family-agnostic. The project's constraint is "never start from a transposase",
which this does not do. TSD is already computed in stage 30 as annotation; the
new part is using it for direction.

TSD ABSENCE IS NOT EVIDENCE OF ANYTHING. IS110/IS1111 leaves no TSD and no TIR
-- measured 0/13 in this project's own Arm B output -- and IS110 is one of the
two dominant length modes here (1279 bp). So a locus with no TSD is simply not
callable by this route and falls through to Layer 3 unchanged. This narrows
Layer 3's domain; it does not replace it.

    92_structural_polarity.py --recon allele_recon/armB_cl0000 \\
        --polarity polarity/all_cl0000 --out structpol
"""
import argparse
import collections
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import revcomp                                    # noqa: E402
from lib.util import log                                         # noqa: E402


def read_alleles(path):
    out, cur = collections.defaultdict(dict), None
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


def find_tsd(long_seq, lcp, ins_len, kmin, kmax):
    """Longest direct repeat flanking the inserted block.

    Two placements are tested because the decomposition can slide by the
    junction ambiguity: the duplicate may sit immediately BEFORE the insert and
    at its END, or at its START and immediately AFTER it. Both are the same
    biological signature read from different sides.
    """
    s, e = lcp, lcp + ins_len
    best = 0
    for k in range(kmax, kmin - 1, -1):
        if s - k < 0 or e + k > len(long_seq):
            continue
        a = long_seq[s - k:s] == long_seq[e - k:e]      # before-start vs end
        b = long_seq[s:s + k] == long_seq[e:e + k]      # start vs after-end
        if a or b:
            best = k
            break
    return best


def find_tir(ins, kmin, kmax):
    """Terminal inverted repeat: insert's head vs revcomp of its tail."""
    best = 0
    for k in range(kmax, kmin - 1, -1):
        if 2 * k > len(ins):
            continue
        if ins[:k] == revcomp(ins[-k:]):
            best = k
            break
    return best


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--polarity", nargs="+", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=15)
    ap.add_argument("--call-min-tsd", type=int, default=4,
                    help="minimum TSD length to make a polarity CALL. A 2-mer "
                         "matches by chance 1/16 of the time and 3-mer 1/64, so "
                         "short repeats are recorded but not trusted.")
    ap.add_argument("--max-tsd", type=int, default=15,
                    help="upper bound on a believable TSD. Real target-site "
                         "duplications are 2-15 bp (IS3 3-4, IS1 9). A 31-60 bp "
                         "direct repeat flanking an insert is NOT a TSD -- it "
                         "means the insert sits inside a larger repeat array or "
                         "is a segmental duplication, and it must not be read as "
                         "polarity. Measured: 9 of 85 calls exceeded 15 bp, "
                         "three of them running past 60.")
    ap.add_argument("--tir-min", type=int, default=8)
    args = ap.parse_args()

    loci, aseq = {}, {}
    for d in args.recon:
        p = os.path.join(d, "loci.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                loci[r["locus_id"]] = r
        aseq.update(read_alleles(os.path.join(d, "alleles.fna")))
    pol = {}
    for d in args.polarity:
        p = os.path.join(d, "polarity.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                pol[r["locus_id"]] = r
    log("%d Layer-2 loci, %d with a Layer-3 record" % (len(loci), len(pol)))

    rows, stat = [], collections.Counter()
    for lid, r in sorted(loci.items()):
        ev = r.get("event_class", "")
        decomposed = ev in ("insertion_target_retained", "replacement",
                            "pure_insertion")
        stat["loci"] += 1
        if not decomposed:
            stat["no_insertion_model"] += 1
            rows.append([lid, ev, 0, 0, ".", "not_decomposed",
                         pol.get(lid, {}).get("polarity_tier", "."),
                         pol.get(lid, {}).get("ancestral_allele", ".")])
            continue
        av = aseq.get(lid, {})
        sid, gid = r.get("shortest_allele"), r.get("longest_allele")
        ss, ls = av.get(sid, ""), av.get(gid, "")
        try:
            lcp = int(r.get("lcp_bp", 0))
            ilen = int(r.get("inserted_len", 0))
        except ValueError:
            continue
        if not ls or ilen <= 0:
            continue
        tsd = find_tsd(ls, lcp, ilen, args.kmin, args.kmax)
        ins = ls[lcp:lcp + ilen]
        tir = find_tir(ins, args.tir_min, 40)
        if tsd > args.max_tsd:
            call, anc = "repeat_array_not_tsd", "."
            stat["oversized_repeat"] += 1
        elif tsd >= args.call_min_tsd:
            call, anc = "derived_is_long", sid
            stat["called_by_tsd"] += 1
        else:
            call, anc = "not_callable", "."
            stat["no_tsd"] += 1
        rows.append([lid, ev, tsd, tir, anc, call,
                     pol.get(lid, {}).get("polarity_tier", "."),
                     pol.get(lid, {}).get("ancestral_allele", ".")])

    with open(args.out + "_calls.tsv", "w") as fh:
        fh.write("locus_id\tevent_class\ttsd_len\ttir_len\tstruct_ancestral_allele\t"
                 "struct_call\tlayer3_tier\tlayer3_ancestral_allele\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    # ---- concordance with Layer 3 P1 ------------------------------------
    agree = dis = 0
    for r in rows:
        lid, _, tsd, _, anc, call, tier, l3 = r
        if call != "derived_is_long" or tier != "P1" or l3 in (".", ""):
            continue
        if anc == l3:
            agree += 1
        else:
            dis += 1

    n = stat["loci"]
    called = stat["called_by_tsd"]
    log("=" * 78)
    log("STRUCTURE-ONLY POLARITY  (TSD >= %d bp, no panel, no tree)"
        % args.call_min_tsd)
    log("")
    log("  Layer-2 loci examined        %5d" % n)
    log("  no insertion model           %5d" % stat["no_insertion_model"])
    log("  decomposed but no TSD        %5d" % stat["no_tsd"])
    log("  oversized repeat, NOT a TSD  %5d" % stat["oversized_repeat"])
    log("  CALLED by TSD                %5d   %5.1f%% of all loci"
        % (called, 100.0 * called / n if n else 0))
    log("")
    tsds = [r[2] for r in rows if r[2] >= args.call_min_tsd]
    if tsds:
        h = collections.Counter(tsds)
        log("  TSD length distribution: " + " ".join(
            "%dbp:%d" % (k, h[k]) for k in sorted(h)))
        log("  median %d bp" % statistics.median(tsds))
    log("")
    log("  CONCORDANCE with Layer 3 P1 (the expensive path):")
    log("    agree     %4d" % agree)
    log("    disagree  %4d" % dis)
    if agree + dis:
        log("    %.1f%% concordant on %d overlapping loci"
            % (100.0 * agree / (agree + dis), agree + dis))
    log("")
    log("  TSD absence is not evidence: IS110 leaves none (0/13 measured here)")
    log("  and is a dominant length mode, so those loci stay with Layer 3.")
    log("=" * 78)
    log("wrote %s_calls.tsv" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
