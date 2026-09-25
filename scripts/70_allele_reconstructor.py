#!/usr/bin/env python3
"""Layer 2 -- reconstruct what the sequence at a locus actually CHANGED FROM.

STATUS: LOCKED 2026-08-28 (assembly tier). place_locus(), exact decomposition,
tolerant decomposition and placement-ambiguity reporting are all locked. Any
change must clear scripts/73_layer2_regression.py, which freezes the 251-event
Arm B benchmark over clusters cl0000/cl0001/cl0002.

WHAT THIS MODULE CLAIMS:
    Given a candidate locus and a cluster of homologous assemblies, reconstruct
    the alternative alleles between homologous anchors, preserve zero-length and
    retained target intervals, and decompose high-homology long/short allele
    pairs into insertion / replacement representations with explicit placement
    and junction ambiguity.

WHAT IT DOES NOT CLAIM:
    A nucleotide-resolved biochemical insertion junction. Arm B's own
    junction_span has median 35 bp (p75 282), so no assembly-only comparison can
    certify +/-2 bp. That certification requires read-level gold and has not
    been done. Do not restate the coarse concordance figures as nucleotide
    precision.

    Polarity is also NOT decided here. Alleles are labelled A, B, C... by
    carrier count; which one is ancestral is Layer 3's problem.

This module answers one question and refuses every other:

    L_ANCHOR | ??? | R_ANCHOR      what is the variable interval, in each genome?

It does NOT decide whether the interval is an IS, a prophage, an Rhs toxin, a
segmental duplication or rRNA. It does not score mobility. Annotation is Layer
3 and runs strictly afterwards, on the reconstructed intervals.

WHY ANCHORS, NEVER A POS. An insertion event is a transformation between two
alleles in a shared frame. That frame is a pair of homologous anchors, not a
coordinate: bubble coordinates live in one genome's space, and the genome that
supplies them is usually a carrier of the element, so coordinates read off it
describe the element's own extent rather than the target site.

WHY THE INTERVAL IS NOT "EMPTY". The pre-event allele is frequently NOT the
concatenation of the two flanks -- it carries a short interval that the
insertion replaced:

    filled : [L_ANCHOR][        INSERTION        ][R_ANCHOR]
    other  : [L_ANCHOR][GTC][R_ANCHOR]

Reconstructing the other allele as left+right silently destroys that GTC, which
is precisely the target-site grammar a later RNA-guide analysis needs. So the
interval between the anchors is extracted verbatim, including when it is 0 bp,
and the transformation recorded is `GTC -> <1310 bp>`, never `empty -> element`.

WHY NO ALLELE IS CALLED "EMPTY". Polarity is not decidable from alleles alone:
`L-X-R` vs `L-R` cannot distinguish insertion from deletion. Alleles are
therefore labelled A, B, C... by descending frequency, and the shortest is
reported as `shortest_allele` with an explicitly frequency-based
`pre_event_confidence`. It is promoted to a directional call ONLY when
--outgroup supplies genomes outside the cluster.

event_class is derived from SEQUENCE TRANSFORMATION ONLY, with no biological
assumption:
    pure_insertion        shortest allele is 0 bp
    replacement           shortest allele is 1..--max-target bp
    length_polymorphism   both alleles long; no clear host/insert relation
    monomorphic           one allele only; nothing happened here
    unresolved            anchors did not place consistently

    70_allele_reconstructor.py --manifest cl0000.manifest \\
        --loci loci.tsv --outdir alleles/cl0000 --anchor 500
"""
import argparse
import csv
import hashlib
import os
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, revcomp, write_fasta                # noqa: E402
from lib import paf as pafmod                                    # noqa: E402
from lib.util import log                                         # noqa: E402
from Bio.Align import PairwiseAligner                            # noqa: E402


def sample_name(path):
    b = os.path.basename(path)
    for suf in (".gz", ".fna", ".fa", ".fasta"):
        if b.endswith(suf):
            b = b[: -len(suf)]
    return b


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()[:12] if s else "EMPTY"


def map_anchors(query_fa, target_fa, out_paf, threads, minimap2):
    """Anchors are ~500 bp and near-identical; asm20 with secondaries kept."""
    cmd = [minimap2, "-c", "-x", "asm20", "-N", "20", "-p", "0.5",
           "--secondary=yes", "-t", str(threads), target_fa, query_fa]
    with open(out_paf, "w") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)
    return out_paf


def best_hits(paf_path, min_ident, min_cov):
    """query name -> list of placements, best first."""
    hits = defaultdict(list)
    for r in pafmod.parse(paf_path):
        if r.identity < min_ident:
            continue
        if r.qspan < min_cov * r.qlen:
            continue
        hits[r.qname].append(r)
    for k in hits:
        hits[k].sort(key=lambda r: (-r.nmatch, -r.identity))
    return hits


def decompose(short, long_):
    """Explain the long allele as short-with-something-inserted.

    The anchors are placed at Arm A's called element edge, which is not the
    true junction, so the variable interval usually retains real target
    sequence on both sides. Measured example: a 32 bp interval where
    short == long[:28] + long[-4:] -- the element landed 28 bp in, and the
    remaining 4 bp sit to its right. Reporting only "32 bp interval" would
    hide that; reporting only "empty site" would destroy it.

    Returns (lcp, lcs, inserted_seq, ambiguity) where
      lcp/lcs    matching prefix / suffix shared by both alleles
      inserted   long[lcp : len(long)-lcs], the sequence that is genuinely new
      ambiguity  lcp + lcs - len(short). > 0 means the junction cannot be
                 placed to the base: that many bp are consistent with either
                 side. This is the microhomology / TSD-like signature, and it
                 is REPORTED, never used to accept or reject the event.
    """
    n = min(len(short), len(long_))
    p = 0
    while p < n and short[p] == long_[p]:
        p += 1
    s = 0
    while s < n - p and short[len(short) - 1 - s] == long_[len(long_) - 1 - s]:
        s += 1
    return p, s, long_[p:len(long_) - s], p + s - len(short)


def _slide(long_, l0, l1):
    """How far can the inserted block slide and still give the same alleles?

    If the base just before the insert equals its last base, the insert can be
    shifted one left with no change to either allele; symmetrically on the
    right. The total is the microhomology at the junction -- REPORTED, never
    resolved, because picking one placement would fake base-level precision.
    """
    left = 0
    while l0 - 1 - left >= 0 and long_[l0 - 1 - left] == long_[l1 - 1 - left]:
        left += 1
    right = 0
    while l1 + right < len(long_) and long_[l0 + right] == long_[l1 + right]:
        right += 1
    return left + right


def tolerant_decompose(short, long_, min_cov, min_ident, min_insert, dominance):
    """Explain long as short-with-an-insert, tolerating SNPs and small indels.

    Exact prefix/suffix matching breaks on the first mismatch, which at
    strain-level divergence (~1 SNP/100 bp) truncates the match after ~80 bp
    and misclassifies most real insertions as length_polymorphism. Measured:
    median (lcp+lcs)/len(short) was 1.00 for cleanly decomposed events but 0.26
    for the failures. So the flanking homology is found by ALIGNMENT instead.

    A strict gate keeps this from manufacturing insertions: the short allele
    must be almost entirely explained (coverage), the explanation must be
    accurate (identity), and there must be ONE dominant extra block rather than
    a scatter of indels (dominance). Anything failing the gate stays
    length_polymorphism / complex -- raising the scorable fraction is not worth
    a single fake insertion.

    Returns a dict, or None when the alignment itself is unusable.
    """
    if not short or not long_ or len(long_) <= len(short):
        return None
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-12, extend_gap_score=-0.4)
    try:
        aln = aligner.align(short, long_)[0]
    except (ValueError, OverflowError, MemoryError):
        return None
    sb, lb = aln.aligned
    if len(sb) == 0:
        return None

    inserts, dels, cols, matches = [], [], 0, 0
    prev_s = prev_l = None
    for (s0, s1), (l0, l1) in zip(sb, lb):
        if prev_s is not None:
            sg, lg = s0 - prev_s, l0 - prev_l
            if sg == 0 and lg > 0:
                inserts.append((lg, prev_s, prev_l, l0))
            elif lg == 0 and sg > 0:
                dels.append(sg)
            elif sg > 0 and lg > 0:
                inserts.append((lg, prev_s, prev_l, l0))
                dels.append(sg)
        seg_s, seg_l = short[s0:s1], long_[l0:l1]
        matches += sum(1 for a, b in zip(seg_s, seg_l) if a == b)
        cols += s1 - s0
        prev_s, prev_l = s1, l1

    if cols == 0:
        return None
    cov = cols / len(short)
    ident = matches / cols
    if not inserts:
        return None
    inserts.sort(key=lambda x: -x[0])
    big, sp, bl0, bl1 = inserts[0]
    others = [x[0] for x in inserts[1:]] + dels
    second = max(others) if others else 0
    n_small_indel = sum(1 for x in others if x > 0)
    dom = big / max(1, second)

    ok = (cov >= min_cov and ident >= min_ident and big >= min_insert
          and (second == 0 or dom >= dominance))
    return {"ok": ok, "coverage": cov, "identity": ident,
            "largest_insert": big, "second_indel": second,
            "n_snp": cols - matches, "n_small_indel": n_small_indel,
            "edit_distance": (cols - matches) + big + sum(others),
            "lcp": sp, "lcs": len(short) - sp,
            "target_lost": sum(dels),
            "insert_seq": long_[bl0:bl1],
            "ambiguity": _slide(long_, bl0, bl1),
            "quality": cov * ident * min(1.0, dom / dominance) if second else cov * ident}


def place_locus(lh, rh, max_interval, prior=None, prior_window=20000,
                margin_ratio=0.05):
    """Score every valid L/R anchor pair and return the best, with a margin.

    The old version maximised raw nmatch subject only to a 50 kb interval cap,
    which on repeat-rich loci let the two anchors land on DIFFERENT copies of a
    repeat. Measured consequence: 3/242 loci placed the wrong interval -- one
    genome handed a 42,367 bp "allele", two others handed the wrong short
    allele entirely. The decomposition was flawless in all three; it was simply
    describing the wrong piece of DNA.

    Deliberately NOT fixed by forcing primary-only alignments (a real insertion
    can legitimately be the secondary hit) nor by constraining the interval to
    the seed's length (an insertion changes that length by design). Instead:

      * anchors normalised by their own length, so a long weak hit cannot
        outvote two clean ones;
      * interval capped at the event scope (--max-interval, now 7 kb, not 50);
      * a LOCALISATION PRIOR on the seed genome only -- Layer 1 already said
        roughly where the locus is, so the seed's anchors have no business
        jumping to a repeat copy tens of kb away. This is positional only and
        never touches junction inference, so it is not circular;
      * the runner-up pair is scored too. If the best pair is not clearly
        better, the placement is reported `ambiguous` rather than forced.

    Returns a dict, or None.
    """
    cands = []
    for L in lh:
        for R in rh:
            if L.tname != R.tname or L.strand != R.strand:
                continue
            if L.strand == "+":
                s, e = L.te, R.ts
            else:
                # anchors swap sides on the minus strand
                s, e = R.te, L.ts
            gap = e - s
            if gap < 0 or gap > max_interval:
                continue
            score = (L.nmatch / max(1, L.qlen)) + (R.nmatch / max(1, R.qlen))
            if prior and L.tname == prior[0]:
                d = max(0, max(prior[1] - e, s - prior[2]))
                if d > prior_window:
                    continue          # seed anchors must stay near their locus
                score -= 0.5 * d / max(1, prior_window)
            cands.append((score, L.tname, s, e, L.strand))
    if not cands:
        return None
    cands.sort(key=lambda x: -x[0])
    best = cands[0]
    second = cands[1][0] if len(cands) > 1 else None
    if second is None:
        margin, status = 1.0, "unique"
    else:
        margin = (best[0] - second) / abs(best[0]) if best[0] else 0.0
        status = "ambiguous" if margin < margin_ratio else "resolved"
    return {"contig": best[1], "start": best[2], "end": best[3],
            "strand": best[4], "n_pairs": len(cands), "best_score": best[0],
            "second_score": second if second is not None else float("nan"),
            "margin": margin, "status": status}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                    help="one genome FASTA path per line (the comparison set)")
    ap.add_argument("--loci", required=True,
                    help="TSV: locus_id, seed_genome, contig, start, end")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--anchor", type=int, default=500,
                    help="homologous anchor length each side (default 500)")
    ap.add_argument("--max-interval", type=int, default=7000,
                    help="reject placements whose interval exceeds this. Sized "
                         "to the <5 kb event scope with headroom, NOT to 50 kb: "
                         "a 50 kb cap let anchors land on different repeat "
                         "copies and fabricate a 42 kb allele.")
    ap.add_argument("--prior-window", type=int, default=20000,
                    help="seed-genome anchors must place within this distance "
                         "of the Layer 1 candidate locus")
    ap.add_argument("--margin-ratio", type=float, default=0.05,
                    help="best anchor pair must beat the runner-up by this "
                         "relative score margin, else placement is ambiguous")
    ap.add_argument("--drop-ambiguous", action="store_true",
                    help="discard ambiguous placements instead of flagging")
    ap.add_argument("--max-target", type=int, default=50,
                    help="shortest allele at/below this is a replaced target "
                         "site rather than a second insertion (default 50)")
    ap.add_argument("--min-ident", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.8,
                    help="fraction of the anchor that must align")
    ap.add_argument("--context", default="50,100,200",
                    help="pre-event target context windows to store")
    ap.add_argument("--outgroup", default=None,
                    help="file of genome names outside the cluster; only with "
                         "this can allele polarity be called directionally")
    ap.add_argument("--min-hom-cov", type=float, default=0.85,
                    help="tolerant gate: fraction of the short allele that must "
                         "be explained as flanking homology")
    ap.add_argument("--min-hom-ident", type=float, default=0.95,
                    help="tolerant gate: identity of that homology")
    ap.add_argument("--min-insert", type=int, default=500,
                    help="MINIMUM INSERT LENGTH, enforced on BOTH decomposition "
                         "paths. Set to 500 on 2026-09-03 from the measured "
                         "mobility cliff, not from a guess about element sizes: "
                         "in K. pneumoniae the within-genome multi-copy rate "
                         "(M1) is 1.3%% at 400-599 bp and 71.8%% at 600-799 bp. "
                         "Below the cliff the only support is M2, at a rate "
                         "consistent with short sequences matching by chance. "
                         "Previously this was a TOLERANT-PATH-ONLY gate at 50 "
                         "while the exact path had NO floor at all, which is "
                         "why inserted_len ran down to 1 bp despite Arm B "
                         "gating bubbles at 300. COST: MITEs (100-400 bp, "
                         "genuinely mobile but non-autonomous) are excluded and "
                         "are out of scope -- on structure alone they are not "
                         "separable from REP arrays and short indels, which is "
                         "what the 0%% M1 rate in those bins reflects.")
    ap.add_argument("--dominance", type=float, default=3.0,
                    help="tolerant gate: largest insert must exceed the next "
                         "largest indel by this factor")
    ap.add_argument("--no-tolerant", action="store_true",
                    help="exact LCP/LCS only (the pre-2026-08-28 behaviour)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ctx_windows = [int(x) for x in args.context.split(",") if x.strip()]
    outg = set()
    if args.outgroup and os.path.exists(args.outgroup):
        outg = {l.strip() for l in open(args.outgroup) if l.strip()}

    paths = [l.strip() for l in open(args.manifest) if l.strip()]
    genomes = {sample_name(p): p for p in paths}
    log("%d genomes in the comparison set" % len(genomes))

    loci = []
    with open(args.loci) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            loci.append({"locus_id": r["locus_id"], "seed": r["seed_genome"],
                         "contig": r["contig"], "start": int(r["start"]),
                         "end": int(r["end"])})
    log("%d candidate loci" % len(loci))

    # ---- build the anchor query set, one L and one R per locus -------------
    fas, tmp = {}, tempfile.mkdtemp(prefix="allrec.")
    qpath = os.path.join(tmp, "anchors.fna")
    kept = []
    with open(qpath, "w") as qf:
        for lo in loci:
            gp = genomes.get(lo["seed"])
            if gp is None:
                continue
            fa = fas.setdefault(lo["seed"], Fasta(gp))
            if lo["contig"] not in fa:
                continue
            clen = fa.length(lo["contig"])
            ls, le = lo["start"] - args.anchor, lo["start"]
            rs, re_ = lo["end"], lo["end"] + args.anchor
            if ls < 0 or re_ > clen:
                continue                      # anchor runs off the contig end
            L = fa.fetch(lo["contig"], ls, le)
            R = fa.fetch(lo["contig"], rs, re_)
            if len(L) < args.anchor or len(R) < args.anchor:
                continue
            write_fasta(qf, lo["locus_id"] + "#L", L)
            write_fasta(qf, lo["locus_id"] + "#R", R)
            kept.append(lo)
    log("%d loci have full-length anchors inside the contig" % len(kept))

    # ---- place the anchors in every genome ---------------------------------
    # locus_id -> genome -> allele sequence
    obs = defaultdict(dict)
    plc = defaultdict(dict)      # locus -> genome -> (contig, s, e, strand)
    pinfo = defaultdict(dict)    # locus -> genome -> placement scoring detail
    amb_drop = {}
    for gname, gp in sorted(genomes.items()):
        p = os.path.join(tmp, "%s.paf" % gname)
        try:
            map_anchors(qpath, gp, p, args.threads, args.minimap2)
        except subprocess.CalledProcessError:
            log("  minimap2 failed on %s, skipped" % gname)
            continue
        hits = best_hits(p, args.min_ident, args.min_cov)
        fa = fas.setdefault(gname, Fasta(gp))
        n = 0
        for lo in kept:
            lh = hits.get(lo["locus_id"] + "#L", [])
            rh = hits.get(lo["locus_id"] + "#R", [])
            if not lh or not rh:
                continue
            prior = ((lo["contig"], lo["start"], lo["end"])
                     if gname == lo["seed"] else None)
            pl = place_locus(lh, rh, args.max_interval, prior,
                             args.prior_window, args.margin_ratio)
            if pl is None:
                continue
            if pl["status"] == "ambiguous" and args.drop_ambiguous:
                amb_drop[lo["locus_id"]] = amb_drop.get(lo["locus_id"], 0) + 1
                continue
            contig, s, e, strand = pl["contig"], pl["start"], pl["end"], pl["strand"]
            pinfo[lo["locus_id"]][gname] = pl
            seq = fa.fetch(contig, s, e) if e > s else ""
            if strand == "-":
                seq = revcomp(seq)
            obs[lo["locus_id"]][gname] = seq
            plc[lo["locus_id"]][gname] = (contig, s, e, strand)
            n += 1
        log("  %-22s placed %d/%d loci" % (gname, n, len(kept)))

    # ---- collapse to alleles ----------------------------------------------
    lo_fh = open(os.path.join(args.outdir, "loci.tsv"), "w")
    al_fh = open(os.path.join(args.outdir, "alleles.tsv"), "w")
    fa_fh = open(os.path.join(args.outdir, "alleles.fna"), "w")
    # NOT "target context". These windows bracket the whole ANCHORED INTERVAL,
    # which when loci are seeded with a pad is ~400 bp of ordinary flanking
    # sequence with the junction buried at offset lcp inside it. Feeding these
    # straight into a motif or RNA-guide analysis puts the window ~200 bp away
    # from the actual target site. Layer 4 re-centres on the junction and
    # re-extracts from the carrier genome; use that instead.
    ctx_fh = open(os.path.join(args.outdir,
                               "anchored_interval_context.tsv"), "w")

    lo_fh.write("\t".join([
        "locus_id", "seed_genome", "contig", "start", "end", "anchor_bp",
        "n_genomes_placed", "n_alleles", "shortest_allele", "shortest_len",
        "longest_allele", "longest_len", "delta_len", "event_class",
        "shortest_n_genomes", "pre_event_confidence", "polarity",
        "allele_len_spread", "lcp_bp", "lcs_bp", "inserted_len",
        "junction_ambiguity_bp", "target_bases_lost", "inserted_md5",
        "target_retained_bp", "n_equivalent_decompositions",
        "decomposition_confidence", "pre_event_genome", "pre_event_contig",
        "pre_event_strand", "junction_L_abs", "junction_R_abs",
        "decomposition_method", "decomposition_quality", "homologous_identity",
        "homologous_coverage", "largest_insert_block_len",
        "second_largest_indel_len", "n_snp", "n_small_indel",
        "alignment_edit_distance", "n_long_alleles_tested",
        "inserted_len_spread", "class_stable",
        "n_candidate_anchor_pairs", "best_pair_score", "second_pair_score",
        "placement_score_margin", "placement_confidence", "placement_status"]) + "\n")
    al_fh.write("\t".join([
        "locus_id", "allele_id", "length", "md5", "n_genomes", "genomes",
        "sequence_if_short"]) + "\n")
    ctx_fh.write("\t".join([
        "locus_id", "allele_id", "window_bp", "source_genome",
        "left_context", "target_interval", "right_context"]) + "\n")

    n_out = 0
    for lo in kept:
        per = obs.get(lo["locus_id"], {})
        if len(per) < 2:
            continue
        groups = defaultdict(list)
        for g, s in per.items():
            groups[s].append(g)
        # allele ids by descending carrier count, then ascending length
        order = sorted(groups.items(), key=lambda kv: (-len(kv[1]), len(kv[0])))
        ids = {}
        for i, (seq, gs) in enumerate(order):
            aid = chr(ord("A") + i) if i < 26 else "A%d" % i
            ids[seq] = aid
            al_fh.write("\t".join([
                lo["locus_id"], aid, str(len(seq)), md5(seq), str(len(gs)),
                ",".join(sorted(gs)),
                seq if len(seq) <= 60 else "."]) + "\n")
            write_fasta(fa_fh, "%s|%s|len=%d|n=%d"
                        % (lo["locus_id"], aid, len(seq), len(gs)), seq or "-")

        lens = sorted(len(s) for s in groups)
        short_seq = min(groups, key=len)
        long_seq = max(groups, key=len)
        sl, ll = len(short_seq), len(long_seq)

        lcp, lcs, ins_seq, amb = decompose(short_seq, long_seq)
        method, tq = "exact", float("nan")
        hid = hcov = float("nan")
        big = second = nsnp = nsi = edist = -1
        td = None
        if not args.no_tolerant and lcp + lcs < len(short_seq):
            # exact decomposition failed to explain the short allele; retry
            # with SNP/indel tolerance before giving up on the event
            td = tolerant_decompose(short_seq, long_seq, args.min_hom_cov,
                                    args.min_hom_ident, args.min_insert,
                                    args.dominance)
        if td is not None:
            hid, hcov = td["identity"], td["coverage"]
            big, second = td["largest_insert"], td["second_indel"]
            nsnp, nsi, edist = td["n_snp"], td["n_small_indel"], td["edit_distance"]
            tq = td["quality"]
            if td["ok"]:
                method = "tolerant_alignment"
                lcp, lcs = td["lcp"], td["lcs"]
                ins_seq, amb = td["insert_seq"], td["ambiguity"]
        # amb >= 0 : junction cannot be placed to the base, that many bp read
        #            equally well on either side (microhomology / TSD-like).
        # amb <  0 : the short allele has bases that survive on NEITHER side of
        #            the junction, i.e. target sequence was genuinely lost.
        ambiguity, lost = max(0, amb), max(0, -amb)
        if method == "exact":
            # BUG FIX 2026-09-02. decompose()'s suffix loop is bounded by
            # `s < n - p` with n = len(short), so lcp + lcs <= len(short) and
            # `amb = lcp + lcs - len(short)` can NEVER be positive on the exact
            # path -- the docstring's "> 0 means the junction cannot be placed
            # to the base" branch was unreachable. Measured consequence: all 57
            # TSD-bearing loci came from the tolerant path and 60 exact loci
            # reported ambiguity 0 while carrying real direct repeats
            # (cl0000.B000000: an 11 bp repeat flanks the insert, reported 0).
            # _slide() computes the true microhomology and is what the tolerant
            # path already used; using it on both paths makes the field
            # reproduce the independent k-sweep detector exactly (76/76 calls,
            # identical lengths).
            ambiguity = _slide(long_seq, lcp, len(long_seq) - lcs)
        tol_lost = td["target_lost"] if (td and method == "tolerant_alignment") else None
        if len(groups) == 1:
            ev = "monomorphic"
        elif sl == ll:
            # every allele the same length: the difference is substitutions, so
            # a prefix/suffix decomposition says nothing. min()/max() would
            # also return the SAME sequence here, silently faking a clean call.
            ev = "substitution_only"
            lcp = lcs = ambiguity = lost = 0
            ins_seq = ""
        elif len(ins_seq) < args.min_insert:
            # THE FLOOR, applied on both paths. The exact path previously had
            # none, so a 1 bp difference decomposed into a reportable "insert".
            ev = "below_min_insert"
            lcp = lcs = ambiguity = lost = 0
        elif sl == 0:
            ev = "pure_insertion"
        elif method == "tolerant_alignment":
            ambiguity, lost = amb, tol_lost or 0
            ev = "replacement" if lost else "insertion_target_retained"
        elif lcp + lcs >= sl:
            # every base of the short allele is accounted for on one side of
            # the junction or the other: a clean insertion into retained target
            ev = "insertion_target_retained"
        elif sl <= args.max_target:
            ev = "replacement"
        else:
            ev = "length_polymorphism"

        n_short = len(groups[short_seq])
        conf = n_short / len(per)
        pol = "unpolarized"
        if outg:
            og_short = sum(1 for g in groups[short_seq] if g in outg)
            og_tot = sum(1 for g in per if g in outg)
            if og_tot:
                pol = ("shortest_is_ancestral" if og_short == og_tot
                       else "shortest_is_derived" if og_short == 0
                       else "outgroup_conflict")

        # FRAME INVARIANCE. Decompose the short allele against every distinct
        # long allele, not just the longest. A stable reconstruction must give
        # the same inserted length and the same class whichever carrier frame
        # it is read from; a carrier-frame bug already bit this module once.
        alt_lens, alt_cls = [], set()
        for cand in groups:
            if cand == short_seq or len(cand) <= sl:
                continue
            p_, s_, ins_, a_ = decompose(short_seq, cand)
            m_ = "exact"
            if not args.no_tolerant and p_ + s_ < sl:
                t_ = tolerant_decompose(short_seq, cand, args.min_hom_cov,
                                        args.min_hom_ident, args.min_insert,
                                        args.dominance)
                if t_ and t_["ok"]:
                    ins_, m_ = t_["insert_seq"], "tolerant_alignment"
                    p_, s_ = t_["lcp"], t_["lcs"]
            alt_lens.append(len(ins_))
            alt_cls.add("replacement" if p_ + s_ < sl and m_ == "exact"
                        else "insertion_target_retained")
        n_alt = len(alt_lens)
        len_spread = (max(alt_lens) - min(alt_lens)) if alt_lens else 0
        stable = "yes" if (n_alt <= 1 or (len_spread == 0 and len(alt_cls) <= 1)) \
            else "no"

        # How many distinct decompositions give the SAME pair of alleles? The
        # insert can slide wherever prefix and suffix agree, so the insertion
        # point is only determined to within `ambiguity` bp. Reported, not
        # resolved: picking one silently would fake base-level precision.
        n_equiv = ambiguity + 1 if ev != "substitution_only" else 0
        dconf = 1.0 / n_equiv if n_equiv else 0.0

        # Absolute junction coordinates in the frame of a genome carrying the
        # pre-event (shortest) allele -- the only frame in which a junction is
        # meaningful, per the Arm B lesson about reading junctions off a carrier.
        pe_g = pe_c = pe_st = "."
        jL_abs = jR_abs = -1
        if ev != "substitution_only" and groups[short_seq]:
            # Prefer the SEED genome as the frame when it carries the shortest
            # allele. The seed defines the locus, so reporting junctions in
            # another carrier's coordinates just makes them incomparable to
            # whatever produced the seed. The junction position itself is still
            # derived from sequence, so this fixes the frame, not the answer.
            carriers = sorted(groups[short_seq])
            pe_g = lo["seed"] if lo["seed"] in groups[short_seq] else carriers[0]
            pp = plc.get(lo["locus_id"], {}).get(pe_g)
            if pp:
                pe_c, ps, pe_, pe_st = pp
                if pe_st == "+":
                    jL_abs, jR_abs = ps + lcp, pe_ - lcs
                else:
                    # allele was revcomp'd: prefix in allele space is the RIGHT
                    # end in contig space
                    jL_abs, jR_abs = ps + lcs, pe_ - lcp

        pi = pinfo.get(lo["locus_id"], {})
        fr = pi.get(pe_g) if pe_g != "." else None
        n_pairs = fr["n_pairs"] if fr else -1
        bsc = "%.4f" % fr["best_score"] if fr else "NA"
        ssc = ("%.4f" % fr["second_score"]
               if fr and fr["second_score"] == fr["second_score"] else "NA")
        pmar = "%.4f" % fr["margin"] if fr else "NA"
        # locus-level status is the WORST across contributing genomes: one
        # ambiguous placement is enough to make the allele set untrustworthy
        sts = [v["status"] for v in pi.values()]
        pstatus = ("ambiguous" if "ambiguous" in sts
                   else ("unique" if sts and all(x == "unique" for x in sts)
                         else ("resolved" if sts else "unknown")))
        pconf = "%.3f" % (min((v["margin"] for v in pi.values()), default=0.0))

        lo_fh.write("\t".join(map(str, [
            lo["locus_id"], lo["seed"], lo["contig"], lo["start"], lo["end"],
            args.anchor, len(per), len(groups), ids[short_seq], sl,
            ids[long_seq], ll, ll - sl, ev, n_short, "%.3f" % conf, pol,
            lens[-1] - lens[0], lcp, lcs, len(ins_seq), ambiguity, lost,
            md5(ins_seq), sl, n_equiv, "%.3f" % dconf,
            pe_g, pe_c, pe_st, jL_abs, jR_abs,
            method, "%.3f" % tq if tq == tq else "NA",
            "%.4f" % hid if hid == hid else "NA",
            "%.4f" % hcov if hcov == hcov else "NA",
            big, second, nsnp, nsi, edist,
            n_alt, len_spread, stable,
            n_pairs, bsc, ssc, pmar, pconf, pstatus])) + "\n")

        # pre-event target context: taken from a genome carrying the SHORTEST
        # allele, because the carrier's own flanks may have been altered by the
        # insertion or by later recombination.
        if ev in ("replacement", "pure_insertion",
                  "insertion_target_retained") and groups[short_seq]:
            src = pe_g if pe_g != "." else sorted(groups[short_seq])[0]
            gp = genomes.get(src)
            if gp:
                fa = fas.setdefault(src, Fasta(gp))
                lh = best_hits(os.path.join(tmp, "%s.paf" % src),
                               args.min_ident, args.min_cov)
                pl = plc.get(lo["locus_id"], {}).get(src)
                if pl:
                    contig, s, e, strand = pl
                    clen = fa.length(contig)
                    for w in ctx_windows:
                        lc = fa.fetch(contig, max(0, s - w), s)
                        rc = fa.fetch(contig, e, min(clen, e + w))
                        ti = fa.fetch(contig, s, e) if e > s else ""
                        if strand == "-":
                            lc, rc = revcomp(rc), revcomp(lc)
                            ti = revcomp(ti)
                        ctx_fh.write("\t".join([
                            lo["locus_id"], ids[short_seq], str(w), src,
                            lc, ti or "-", rc]) + "\n")
        n_out += 1

    for fh in (lo_fh, al_fh, fa_fh, ctx_fh):
        fh.close()

    # ---- summary -----------------------------------------------------------
    cls = Counter()
    with open(os.path.join(args.outdir, "loci.tsv")) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        cls[r["event_class"]] += 1
    log("=" * 72)
    log("ALLELE RECONSTRUCTION: %d loci resolved in >=2 genomes" % n_out)
    for k, v in cls.most_common():
        log("  %-22s %5d" % (k, v))
    rec = [r for r in rows
           if r["event_class"] in ("insertion_target_retained", "replacement",
                                   "pure_insertion")]
    if rec:
        il = [int(r["inserted_len"]) for r in rec]
        tl = [int(r["target_retained_bp"]) for r in rec]
        am = [int(r["junction_ambiguity_bp"]) for r in rec]
        lo_ = [int(r["target_bases_lost"]) for r in rec]
        log("")
        log("  inserted length      : median %d bp, range %d-%d"
            % (statistics.median(il), min(il), max(il)))
        log("  target retained      : median %d bp, range %d-%d"
            % (statistics.median(tl), min(tl), max(tl)))
        log("  junction ambiguity   : median %d bp, range %d-%d  "
            "(microhomology / TSD-like; reported, never a filter)"
            % (statistics.median(am), min(am), max(am)))
        log("  target bases lost    : median %d bp, range %d-%d"
            % (statistics.median(lo_), min(lo_), max(lo_)))
        hist = Counter(il)
        log("  inserted-length modes: " + " ".join(
            "%dbp:%d" % (k, v) for k, v in hist.most_common(8)))
    meth = Counter(r["decomposition_method"] for r in rows)
    stab = Counter(r["class_stable"] for r in rows)
    log("")
    log("  decomposition method: " + ", ".join("%s=%d" % kv for kv in meth.most_common()))
    log("  frame-invariant     : " + ", ".join("%s=%d" % kv for kv in stab.most_common()))
    pst = Counter(r["placement_status"] for r in rows)
    log("  placement status    : " + ", ".join("%s=%d" % kv for kv in pst.most_common()))
    log("=" * 72)
    log("wrote %s/{loci,alleles,anchored_interval_context}.tsv + alleles.fna"
        % args.outdir)


if __name__ == "__main__":
    main()
