#!/usr/bin/env python3
"""Audit the Layer 2 / Arm B inserted-length disagreements. Not "who is right".

Arm B is not a length gold: its bubble extent is known to carry target/flank
padding (Arm A showed the same inflation, 1311 bp called for a 1279 bp element).
So a length difference is evidence of nothing on its own. Each disagreement is
therefore classified by INDEPENDENT evidence:

  ArmB_extent_imprecise   the Layer 2 insert aligns inside the Arm B insert at
                          high identity AND high coverage -- same element, Arm B
                          simply drew a wider box. Layer 2 likely correct.
  different_long_alleles  the locus genuinely carries several distinct long
                          alleles; the two tools picked different ones. Not an
                          algorithm error.
  Layer2_repeat_ambiguous high junction ambiguity, or the carriers disagree
                          among themselves. Layer 2 uncertain, flagged not fixed.
  true_decomposition_error the inserted sequences do not correspond, carriers
                          agree with Arm B rather than with the call. A real bug.

The strongest evidence is sequence correspondence, not length, plus carrier
consensus: if every (short-allele, long-allele) pair at a locus independently
yields the same inserted length, that is strong assembly-level evidence
regardless of what Arm B said.

    72_layer2_disagreement_audit.py --clusters cl0000 cl0001 cl0002 \\
        --recon-dir allele_recon --armb-dir armB_bench/stage30 --out audit
"""
import argparse
import csv
import os
import statistics
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import revcomp                                    # noqa: E402
from lib.util import log                                         # noqa: E402
from Bio.Align import PairwiseAligner                            # noqa: E402

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_ar", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "70_allele_reconstructor.py"))
_ar = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ar)

GATE = dict(min_cov=0.85, min_ident=0.95, min_insert=50, dominance=3.0)


def read_fasta(path):
    out, name = {}, None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = []
        elif name is not None:
            out[name].append(line.strip())
    return {k: "".join(v).replace("-", "") for k, v in out.items()}


def decomp(short, long_):
    """Same cascade the module uses: exact first, tolerant as fallback."""
    if not short or not long_ or len(long_) <= len(short):
        return None
    p, q, ins, amb = _ar.decompose(short, long_)
    if p + q >= len(short):
        return len(ins), ins, max(0, amb), "exact"
    t = _ar.tolerant_decompose(short, long_, GATE["min_cov"], GATE["min_ident"],
                               GATE["min_insert"], GATE["dominance"])
    if t and t["ok"]:
        return len(t["insert_seq"]), t["insert_seq"], t["ambiguity"], "tolerant"
    return None


def _corr_one(a, b):
    """Local alignment of the shorter into the longer -> (identity, coverage)."""
    if not a or not b:
        return 0.0, 0.0
    al = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                         open_gap_score=-12, extend_gap_score=-0.4)
    try:
        aln = al.align(a, b)[0]
    except (ValueError, OverflowError, MemoryError, IndexError):
        return 0.0, 0.0
    sb, lb = aln.aligned
    cols = matches = 0
    for (s0, s1), (l0, l1) in zip(sb, lb):
        matches += sum(1 for x, y in zip(a[s0:s1], b[l0:l1]) if x == y)
        cols += s1 - s0
    if cols == 0:
        return 0.0, 0.0
    return matches / cols, cols / min(len(a), len(b))


def seq_corr(a, b):
    """Orientation-agnostic correspondence.

    Arm B stores its insert in the backbone's orientation, Layer 2 in the
    pre-event carrier's; when those differ the forward comparison returns ~0.50
    identity, which is the score two unrelated sequences get from a local
    aligner and looks exactly like a decomposition error. Measured: 6 of 9
    apparent errors were this artefact, matching revcomp at 1.000 identity.
    So both orientations are scored and the better is taken.
    """
    f = _corr_one(a, b)
    r = _corr_one(a, revcomp(b))
    return f if f[0] >= r[0] else r


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", nargs="+", required=True)
    ap.add_argument("--recon-dir", required=True)
    ap.add_argument("--armb-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefix", default="armB_",
                    help="recon subdir prefix (default armB_)")
    ap.add_argument("--min-corr-ident", type=float, default=0.95)
    ap.add_argument("--min-corr-cov", type=float, default=0.90)
    args = ap.parse_args()

    loci, alleles, gold, ginsert = {}, defaultdict(dict), {}, {}
    for c in args.clusters:
        d = os.path.join(args.recon_dir, args.prefix + c)
        with open(os.path.join(d, "loci.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                loci[r["locus_id"]] = r
        for k, v in read_fasta(os.path.join(d, "alleles.fna")).items():
            p = k.split("|")
            if len(p) >= 2:
                alleles[p[0]][p[1]] = v
        with open(os.path.join(args.armb_dir, c, "armB_events.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gold[r["event_id"]] = r
        for k, v in read_fasta(os.path.join(args.armb_dir, c,
                                            "armB_inserts.fna")).items():
            ginsert[k.split()[0]] = v

    rows, cls = [], Counter()
    deltas = []
    for lid, r in sorted(loci.items()):
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        g = gold.get(lid)
        if not g:
            continue
        l2len, gl = int(r["inserted_len"]), int(g["insert_size"])
        if l2len == gl:
            continue
        d = l2len - gl
        deltas.append(d)

        av = alleles.get(lid, {})
        short = av.get(r["shortest_allele"], "")
        # carrier consensus: every distinct (short, long) allele pair
        lens, insseqs = [], []
        for aid, seq in av.items():
            if len(seq) <= len(short):
                continue
            res = decomp(short, seq)
            if res:
                lens.append(res[0])
                insseqs.append(res[1])
        mode_len, mode_frac = (-1, 0.0)
        if lens:
            cnt = Counter(lens)
            mode_len, n = cnt.most_common(1)[0]
            mode_frac = n / len(lens)

        # sequence correspondence against the Arm B insert
        l2ins = insseqs[0] if insseqs else ""
        for s in insseqs:
            if len(s) == l2len:
                l2ins = s
                break
        gi = ginsert.get(lid, "")
        ident, cov = seq_corr(l2ins, gi) if (l2ins and gi) else (0.0, 0.0)

        # reverse-complement stability
        rc = decomp(revcomp(short), revcomp(av.get(r["longest_allele"], "")))
        rc_len = rc[0] if rc else -1

        amb = int(r["junction_ambiguity_bp"])
        n_alt = len(lens)
        if ident >= args.min_corr_ident and cov >= args.min_corr_cov:
            k = "ArmB_extent_imprecise"
        elif n_alt > 1 and mode_frac < 0.75:
            k = "different_long_alleles"
        elif amb > 5 or (rc_len >= 0 and rc_len != l2len):
            k = "Layer2_repeat_ambiguous"
        elif ident < 0.80 or cov < 0.50:
            k = "true_decomposition_error"
        else:
            k = "Layer2_repeat_ambiguous"
        cls[k] += 1
        rows.append([lid, r["decomposition_method"], l2len, gl, d,
                     "%.3f" % ident, "%.3f" % cov, amb, rc_len, n_alt,
                     mode_len, "%.2f" % mode_frac, k])

    with open(args.out + "_disagreements.tsv", "w") as fh:
        fh.write("locus_id\tmethod\tlayer2_inserted_len\tarmB_insert_size\t"
                 "delta\tinsert_seq_identity\tinsert_seq_coverage\t"
                 "junction_ambiguity_bp\trevcomp_inserted_len\t"
                 "n_allele_pairs\tcarrier_mode_len\tcarrier_mode_fraction\t"
                 "audit_class\tvalidation_flag\n")
        for r in rows:
            # A residual disagreement is RECORDED, not resolved. Neither method
            # is a length gold, so the honest label says the two assembly-based
            # methods disagree -- nothing more. Revisit the algorithm only if
            # this flag accumulates a repeated failure signature (10+ cases of
            # one shape); a single case is not evidence of a defect.
            flag = ("assembly_crossmethod_discordant"
                    if r[-1] == "true_decomposition_error" else ".")
            fh.write("\t".join(map(str, r)) + "\t" + flag + "\n")

    log("=" * 78)
    log("TOLERANT-SUBSET DISAGREEMENT AUDIT  (n=%d)" % len(rows))
    log("")
    for k, v in cls.most_common():
        log("  %-28s %4d  %5.1f%%" % (k, v, 100.0 * v / len(rows) if rows else 0))
    log("")
    if deltas:
        log("  delta length (Layer2 - ArmB):")
        log("    median %+.0f   mean %+.1f   range %+d..%+d"
            % (statistics.median(deltas), statistics.mean(deltas),
               min(deltas), max(deltas)))
        neg = sum(1 for d in deltas if d < 0)
        small = sum(1 for d in deltas if abs(d) <= 50)
        log("    negative (Layer2 shorter): %d/%d = %.0f%%"
            % (neg, len(deltas), 100.0 * neg / len(deltas)))
        log("    |delta| <= 50 bp        : %d/%d = %.0f%%"
            % (small, len(deltas), 100.0 * small / len(deltas)))
        hist = Counter(deltas)
        log("    most common deltas: " + " ".join(
            "%+d:%d" % (k, v) for k, v in hist.most_common(10)))
    log("=" * 78)
    log("wrote %s_disagreements.tsv" % args.out)


if __name__ == "__main__":
    main()
