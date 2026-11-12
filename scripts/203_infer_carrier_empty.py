#!/usr/bin/env python3
"""Infer, per carrier occurrence, the site with that insertion removed.

WHAT THIS PRODUCES. inferred_empty_site = the carrier's OWN upstream joined
directly to the carrier's OWN downstream, across the located insert interval.
Every base comes from that carrier assembly, so the carrier's SNPs and indels
are retained. The observed empty allele is used ONLY to identify the homologous
locus and to judge junction structure -- never to substitute bases, and never
to fill a core the carrier cannot supply. It is an INFERENCE in that carrier's
background, not a claim about an ancestral sequence. Nothing is overwritten:
catalogue_v3, database_v3 and insertions_all/ are read-only here.

THE UNIT IS AN OCCURRENCE:
    event_id + carrier_accession + carrier_contig + insert_start + insert_end
Events are deduplicated on sequence, so one insert may occur in many carriers
and at several loci in one carrier. Each occurrence is resolved and written
separately; an inferred empty site is NEVER shared between occurrences.

COORDINATES. 0-based, end-exclusive: carrier_contig[s:e] is the insert, so
e - s == inserted_len. The empty-side junction_L_abs/junction_R_abs are also
0-based slice bounds, verified in that frame. The two frames use the SAME
FORMAT but are NEVER interchangeable as VALUES: s/e come from locating the
insert in the carrier itself. A 1-based inclusive rendering of [s,e) is
[s+1, e].

WHY NO EXTRA TRIM -- AND WHAT THAT CLAIM RESTS ON.
An earlier version of this file argued the insert interval contains one TSD
copy because (len_long - len_short) - inserted_len was 0 for 77% of events.
THAT ARGUMENT IS WITHDRAWN. decompose() defines inserted = long[lcp:-lcs], so
that difference IS decompose's own ambiguity lcp+lcs-len(short), which is 0 for
every clean decomposition (64,809/64,809 on the exact path). It was an identity
the code guarantees, not evidence. Separately, junction_overlap_bp is NOT a
verified TSD: it is _slide(), "the microhomology at the junction -- REPORTED,
never resolved", i.e. an alignment-derived overlap candidate.
What does hold: sliding is sequence-invariant by construction, so excising
inserted_len bases at ANY valid slide position reproduces the empty allele.
Hence the DEFAULT is upstream + downstream with no extra trimming, and a second
TSD copy is never removed. But slide-invariance is a property of the ALLELE; it
does not establish that the located CARRIER position is a valid slide position.
That is what the junction checks below test, and why reconstruction_status is
reported separately from localization_status.

TWO STATUSES, BECAUSE THEY ARE DIFFERENT QUESTIONS.
  localization_status    is this carrier copy the one corresponding to the
                         event's locus?
  reconstruction_status  can a clean empty sequence be built at that copy?
A sequence can be constructible at a copy whose correspondence is unresolved,
and vice versa. Only occurrences resolved on BOTH reach the official FASTA.

IDENTITY, COVERAGE, MARGIN -- EXACT FORMULAS.
  aligned columns  the ungapped overlap of the anchor and the carrier flank at
                   the best shift in [-MAX_SHIFT, +MAX_SHIFT]
  identity   = matches / aligned_columns          (0 if aligned_columns == 0)
  coverage   = aligned_columns / anchor_bp
  margin     = min(up_id,dn_id) of the best candidate
               minus min(up_id,dn_id) of the runner-up
A nonzero best shift means the carrier boundary sits that many bp from where
the empty-side anchor says it should: positive => insert bases would be left
behind, negative => host bases would be deleted. High overall identity hides
this, which is why the shift is reported per side rather than folded into a
score.

CANDIDATE COMPLETENESS. Copies are enumerated until EXAMINE_CAP. If the cap is
reached the search is marked candidate_search_truncated and the occurrence may
NOT be called uniquely localized, because a better copy may be unexamined.
Several hits are candidates, never confirmed occurrences.
"""
import argparse
import collections
import hashlib
import os
import sqlite3
import sys
import time

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"
OUT = W + "/inferred_empty_v1"
DB = W + "/locations/insertions.sqlite"
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

CAND_MIN = 0.80          # keep as a candidate
TRUST_ID = 0.90          # trusted identity, per side (configurable)
TRUST_COV = 0.90         # trusted coverage, per side (configurable)
AMBIG_MARGIN = 0.02
MAX_SHIFT = 10           # junction displacement searched, each direction
EXAMINE_CAP = 500        # copies examined before declaring truncation


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    out, k, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if k is not None:
                    out[k] = "".join(buf).upper()
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if k is not None:
        out[k] = "".join(buf).upper()
    return out


def load_genome(path, tries=4):
    for t in range(tries):
        try:
            return read_fasta(path)
        except (IOError, OSError):
            if t == tries - 1:
                return None
            time.sleep(0.5 * (t + 1))
    return None


def load_map(path):
    m = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                m[f[0]] = f[1]
    return m


def cmp_ungapped(a, b):
    """(identity, aligned_columns, matches) over the right-anchored overlap."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0, 0, 0
    m = sum(1 for x, y in zip(a[-n:], b[-n:]) if x == y)
    return m / float(n), n, m


def flank_at(cseq, pos, side, shift, A):
    """Carrier flank abutting a boundary placed `shift` bp from `pos`.

    The boundary is moved ON THE CONTIG rather than trimming a pre-cut flank,
    so the returned sequence is always A bp of real carrier context ending
    (up) or starting (dn) exactly at the hypothesised boundary.
    """
    q = pos + shift
    if side == "up":
        return cseq[max(0, q - A):q]
    return cseq[q:q + A]


def best_pairing(cseq, s, e, up_a, dn_a, A, max_shift=MAX_SHIFT):
    """Best (orientation, (up_id,up_aln,up_shift), (dn_id,dn_aln,dn_shift)).

    Orientation is a property of the WHOLE locus, so "same" and "flip" are
    evaluated as units: under flip the empty-side upstream anchor corresponds
    to the reverse complement of the carrier's RIGHT flank, and the downstream
    anchor to the reverse complement of the carrier's LEFT flank. Scoring each
    side independently against whichever carrier flank fits best would let a
    candidate pass with left matched to left and right matched to left, which
    is not a consistent locus.
    """
    def scan(anchor, pos, side, rc_it):
        best = (0.0, 0, 0)
        for sh in range(-max_shift, max_shift + 1):
            seg = flank_at(cseq, pos, side, sh, A)
            if rc_it:
                seg = revcomp(seg)
            idn, aln, _ = cmp_ungapped(anchor, seg)
            if aln and (idn, aln) > (best[0], best[1]):
                best = (idn, aln, sh)
        return best

    same = (scan(up_a, s, "up", False), scan(dn_a, e, "dn", False))
    flip = (scan(up_a, e, "dn", True), scan(dn_a, s, "up", True))
    if min(same[0][0], same[1][0]) >= min(flip[0][0], flip[1][0]):
        return "same", same[0], same[1]
    return "flip", flip[0], flip[1]


def occurrences(ctgs, needle, cap=EXAMINE_CAP):
    """(hits, truncated). Hits are (contig, start, strand), 0-based."""
    out = []
    r = revcomp(needle)
    for k, v in ctgs.items():
        for pat, st in ((needle, "+"), (r, "-")):
            i = v.find(pat)
            while i >= 0:
                out.append((k, i, st))
                if len(out) >= cap:
                    return out, True
                i = v.find(pat, i + 1)
    return out, False


COLS = [
    "occurrence_id", "db_id", "species", "event_id",
    "carrier_accession", "carrier_contig", "carrier_contig_len",
    "insert_start", "insert_end", "insert_strand", "inserted_len", "insert_md5",
    "empty_accession", "empty_contig", "empty_junction_L", "empty_junction_R",
    "n_candidates_examined", "candidate_search_truncated", "n_candidates_trusted",
    "localization_status", "orientation", "anchor_bp",
    "up_identity", "up_coverage", "up_aligned", "up_shift",
    "dn_identity", "dn_coverage", "dn_aligned", "dn_shift",
    "margin",
    "reconstruction_status", "boundary_evidence",
    "residual_insert_bp", "host_overdeleted_bp", "tsd_state", "core_recoverable",
    "carrier_up_len", "carrier_dn_len",
    "inferred_empty_len", "observed_empty_len", "len_delta", "len_delta_class",
    "consistency_vs_observed_id",
    "junction_overlap_bp", "decomposition_method", "uncertainty_reason",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=SRC + "/genome_map.tsv")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--flank", type=int, default=2000, help="bp of carrier context kept")
    ap.add_argument("--anchor", type=int, default=500, help="bp used to test homology")
    ap.add_argument("--trust-id", type=float, default=TRUST_ID)
    ap.add_argument("--trust-cov", type=float, default=TRUST_COV)
    ap.add_argument("--cand-min", type=float, default=CAND_MIN)
    ap.add_argument("--examine-cap", type=int, default=EXAMINE_CAP)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out + "/parts", exist_ok=True)

    man = load_map(a.manifest)
    ins = read_fasta(SRC + "/insert.fna")
    obs = read_fasta(SRC + "/empty_site.fna")
    meta = {}
    with open(SRC + "/insertions.tsv") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            meta[f[ix["db_id"]]] = (f[ix["junction_overlap_bp"]],
                                    f[ix["decomposition_method"]])

    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    ev = {}
    for r in con.execute(
            "SELECT db_id, species, event_id, coord_frame, coord_contig,"
            "  junction_L_abs, junction_R_abs, insert_md5"
            " FROM loc WHERE junction_L_abs NOT IN ('.','') GROUP BY db_id"):
        ev[r[0]] = r
    by_car = collections.defaultdict(list)
    for db_id, g in con.execute(
            "SELECT db_id, genome FROM loc WHERE role='carrier_derived'"):
        if db_id in ev:
            by_car[g].append(db_id)

    gs = sorted(by_car)
    mine = [g for i, g in enumerate(gs) if i % a.nshard == a.shard]
    sys.stderr.write("shard %d/%d: %d carriers, %d (event,carrier) pairs\n"
                     % (a.shard, a.nshard, len(mine),
                        sum(len(by_car[g]) for g in mine)))

    # empty-side anchors, one read per empty assembly
    need = collections.defaultdict(list)
    for g in mine:
        for db_id in by_car[g]:
            need[ev[db_id][3]].append(db_id)
    A = a.anchor
    anch = {}
    for n, (gf, ids) in enumerate(need.items(), 1):
        p = man.get(gf)
        ctgs = load_genome(p) if p else None
        for db_id in ids:
            r = ev[db_id]
            seq = ctgs.get(r[4]) if ctgs else None
            if seq is None:
                anch[db_id] = None
                continue
            jl, jr = int(r[5]), int(r[6])
            anch[db_id] = (seq[max(0, jl - A):jl], seq[jr:jr + A])
        if n % 300 == 0:
            sys.stderr.write("  anchors %d/%d\n" % (n, len(need)))
    sys.stderr.write("  anchors ready (%d)\n" % len(anch))

    ftsv = open("%s/parts/occ_%03d.tsv" % (a.out, a.shard), "w")
    ffna = open("%s/parts/inferred_%03d.fna" % (a.out, a.shard), "w")
    fother = open("%s/parts/unresolved_%03d.fna" % (a.out, a.shard), "w")
    fcand = open("%s/parts/candidates_%03d.tsv" % (a.out, a.shard), "w")
    stat = collections.Counter()

    def emit(row, seq=None, official=False):
        ftsv.write("\t".join(str(row.get(c, ".")) for c in COLS) + "\n")
        if seq is not None:
            fh = ffna if official else fother
            fh.write(">%s loc=%s rec=%s up=%d dn=%d\n"
                     % (row["occurrence_id"], row["localization_status"],
                        row["reconstruction_status"], row["carrier_up_len"],
                        row["carrier_dn_len"]))
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    for n, g in enumerate(mine, 1):
        p = man.get(g)
        ctgs = load_genome(p) if p else None
        for db_id in by_car[g]:
            r = ev[db_id]
            ov, dmeth = meta.get(db_id, (".", "."))
            base = {"db_id": db_id, "species": r[1], "event_id": r[2],
                    "carrier_accession": g, "empty_accession": r[3],
                    "empty_contig": r[4], "empty_junction_L": r[5],
                    "empty_junction_R": r[6], "insert_md5": r[7],
                    "anchor_bp": A, "junction_overlap_bp": ov,
                    "decomposition_method": dmeth,
                    "occurrence_id": "%s|%s|.|.|." % (db_id, g),
                    "carrier_up_len": 0, "carrier_dn_len": 0}
            if ctgs is None:
                base.update(localization_status="no_carrier_genome",
                            reconstruction_status="not_attempted",
                            uncertainty_reason="carrier assembly unreadable")
                stat["no_carrier_genome"] += 1
                emit(base)
                continue
            iseq = ins.get(db_id, "")
            an = anch.get(db_id)
            if not iseq:
                base.update(localization_status="no_insert_sequence",
                            reconstruction_status="not_attempted",
                            uncertainty_reason="no insert sequence for event")
                stat["no_insert_sequence"] += 1
                emit(base)
                continue
            if an is None:
                base.update(localization_status="no_empty_anchor",
                            reconstruction_status="not_attempted",
                            uncertainty_reason="empty assembly/contig unavailable")
                stat["no_empty_anchor"] += 1
                emit(base)
                continue
            up_a, dn_a = an
            hits, trunc = occurrences(ctgs, iseq, a.examine_cap)
            base["n_candidates_examined"] = len(hits)
            base["candidate_search_truncated"] = "yes" if trunc else "no"
            if not hits:
                base.update(localization_status="insert_not_in_carrier",
                            reconstruction_status="not_attempted",
                            uncertainty_reason="no exact copy of the insert")
                stat["insert_not_in_carrier"] += 1
                emit(base)
                continue

            scored = []
            for ctg, s, st in hits:
                cs = ctgs[ctg]
                e = s + len(iseq)
                o, u, d = best_pairing(cs, s, e, up_a, dn_a, A)
                scored.append((min(u[0], d[0]), u, d, o, ctg, s, st))
            scored.sort(key=lambda x: -x[0])
            best = scored[0]
            margin = best[0] - scored[1][0] if len(scored) > 1 else 1.0
            _, u, d, orient, ctg, s, st = best
            e = s + len(iseq)
            cs = ctgs[ctg]
            up_cov = u[1] / float(A) if A else 0.0
            dn_cov = d[1] / float(A) if A else 0.0
            n_trusted = sum(1 for x in scored
                            if min(x[1][0], x[2][0]) >= a.trust_id)
            base.update(
                occurrence_id="%s|%s|%s|%d|%d" % (db_id, g, ctg, s, e),
                carrier_contig=ctg, carrier_contig_len=len(cs),
                insert_start=s, insert_end=e, insert_strand=st,
                inserted_len=len(iseq), orientation=orient,
                n_candidates_trusted=n_trusted,
                up_identity="%.4f" % u[0], up_coverage="%.4f" % up_cov,
                up_aligned=u[1], up_shift=u[2],
                dn_identity="%.4f" % d[0], dn_coverage="%.4f" % dn_cov,
                dn_aligned=d[1], dn_shift=d[2],
                margin="%.4f" % margin)

            trusted = (u[0] >= a.trust_id and d[0] >= a.trust_id
                       and up_cov >= a.trust_cov and dn_cov >= a.trust_cov)
            candidate = (u[0] >= a.cand_min and d[0] >= a.cand_min)
            reasons = []
            if trunc:
                reasons.append("candidate_search_truncated")
            if not candidate:
                # ONE-SIDED HOMOLOGY IS ITS OWN CATEGORY, NOT GENERIC FAILURE.
                # In the sample 604 of 709 unsupported pairs had exactly one
                # anchor >= the candidate floor while the other was near
                # background. That is the signature of a locus sharing one
                # flank only -- a nested or tandem arrangement, or a structural
                # difference abutting the site -- and it is a different claim
                # from "no homologous locus found", so the two are not merged.
                hi = max(u[0], d[0])
                lo_ = min(u[0], d[0])
                if hi >= a.trust_id and lo_ < a.cand_min:
                    loc = "one_sided_homology"
                    reasons.append("only one anchor supports the locus "
                                   "(%.3f vs %.3f)" % (hi, lo_))
                else:
                    loc = "no_homologous_locus"
                    reasons.append("both anchors below %.2f" % a.cand_min)
            elif not trusted:
                loc = "low_confidence"
                reasons.append("anchor identity/coverage below trust threshold")
            elif trunc:
                loc = "ambiguous_copy"
                reasons.append("cap reached; a better copy may be unexamined")
            elif n_trusted > 1 and margin < AMBIG_MARGIN:
                loc = "ambiguous_copy"
                reasons.append("%d trusted copies within margin %.3f"
                               % (n_trusted, margin))
            else:
                loc = "resolved"
            base["localization_status"] = loc
            stat[loc] += 1

            # ---- build the candidate sequence (default: no extra trim)
            F = a.flank
            up = cs[max(0, s - F):s]
            dn = cs[e:e + F]
            inferred = up + dn
            base["carrier_up_len"] = len(up)
            base["carrier_dn_len"] = len(dn)
            observed = obs.get(db_id, "")
            lo = min(len(inferred), len(observed))
            cid = cmp_ungapped(inferred[:lo], observed[:lo])[0] if lo else 0.0
            base["inferred_empty_len"] = len(inferred)
            base["observed_empty_len"] = len(observed)
            base["len_delta"] = len(inferred) - len(observed)
            base["consistency_vs_observed_id"] = "%.4f" % cid
            # fixed-length windows make len_delta uninformative about the splice
            base["len_delta_class"] = ("fixed_window_not_comparable"
                                       if len(up) == F and len(dn) == F
                                       else "unexplained_length_difference")

            # ---- junction-level evidence
            resid = max(0, u[2]) + max(0, -d[2])
            overdel = max(0, -u[2]) + max(0, d[2])
            core_ok = (len(up) >= min(F, s) and len(dn) >= min(F, len(cs) - e))
            base["residual_insert_bp"] = resid
            base["host_overdeleted_bp"] = overdel
            base["core_recoverable"] = "yes" if core_ok else "no"
            try:
                ovn = int(ov)
            except ValueError:
                ovn = 0
            base["tsd_state"] = ("slide_span_%d_not_verified_TSD" % ovn
                                 if ovn else "no_microhomology_reported")
            if u[2] == 0 and d[2] == 0:
                base["boundary_evidence"] = "both_anchors_abut_insert_boundary"
                rec = "resolved"
            else:
                base["boundary_evidence"] = ("anchor_shift_up=%d_dn=%d"
                                             % (u[2], d[2]))
                rec = "boundary_ambiguous"
                reasons.append("junction displaced (up %+d, dn %+d)" % (u[2], d[2]))
            if not core_ok:
                rec = "core_not_recoverable"
                reasons.append("contig edge truncates a flank")
            if loc in ("one_sided_homology", "no_homologous_locus",
                       "low_confidence"):
                rec = rec if rec != "resolved" else "built_but_locus_unresolved"
            base["reconstruction_status"] = rec
            base["uncertainty_reason"] = "; ".join(reasons) if reasons else "."
            stat["rec_" + rec] += 1

            if loc == "ambiguous_copy" or n_trusted > 1:
                for sc in scored[:20]:
                    fcand.write("%s\t%s\t%s\t%d\t%s\t%.4f\t%.4f\t%s\t%d\t%d\n"
                                % (db_id, g, sc[4], sc[5], sc[6],
                                   sc[1][0], sc[2][0], sc[3], sc[1][2], sc[2][2]))
            official = (loc == "resolved" and rec == "resolved")
            stat["official" if official else "withheld"] += 1
            emit(base, inferred, official)
        if n % 200 == 0:
            sys.stderr.write("  %d/%d carriers\n" % (n, len(mine)))

    for fh in (ftsv, ffna, fother, fcand):
        fh.close()
    print("shard %d: %s" % (a.shard, dict(stat)))


if __name__ == "__main__":
    main()
