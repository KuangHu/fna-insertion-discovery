#!/usr/bin/env python3
"""Stage 2 of the v1 record-format export: one main record per occurrence.

Three passes, so each assembly is opened once:
  --pass carrier    grouped by carrier accession -> insertion, upstream,
                    downstream, reconstructed_empty
  --pass reference  grouped by the reference THIS ROW used -> reference_empty
  --pass merge      joins them into occurrences.jsonl + the view tables

DECISIONS THIS EXPORT MAKES, RECORDED RATHER THAN IMPLIED
---------------------------------------------------------
chosen_joint_shift_bp IS ALWAYS 0. v2's `chosen_shift` is a REPORTED finding
("a nonzero shift was verified sequence-equivalent"), not an APPLIED one: the
published bases were cut at insert_start/insert_end, and a sampled replay
confirms it. Writing the nonzero value into chosen_joint_shift_bp would make
final_start = original + k disagree with the exported upstream flank by k
bases, so the nonzero value is preserved under its own key instead.

THE REFERENCE WINDOW IS DEFINED HERE, NOT INHERITED. v2 never stored a
reference window; it stored a junction. The window is therefore the carrier's
own flank lengths mirrored about the reference junction and clipped to the
reference contig, which makes reference_empty and reconstructed_empty directly
comparable. junction_offset0 = junction0 - window_start0 always, which is why
the window bounds and the junction are separate fields.

FLANK LENGTHS ARE READ PER ROW. They are nominally 2000 but clipped at contig
ends; 55 rows are short. Nothing here assumes 2000.
"""
import argparse
import collections
import hashlib
import json
import os
import sys
import zlib

W = "/global/scratch/users/kh36969/fna_ins_discovery"
SRC = W + "/insertions_all"
V2 = W + "/inferred_empty_v2"
OUT = W + "/insertion_db_v1"
SCHEMA = "1.0"
RECIPE = "concat_carrier_flanks_v2"
INTERVAL_DEF = "catalogue_excised_allele_may_include_extra_TSD"
ROLES = ("insertion", "upstream", "downstream", "reference_empty",
         "reconstructed_empty")


def sid(acc, rec):
    return "SRC|%s|%s" % (acc, rec)


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def fa(p):
    o, k, b = {}, None, []
    with open(p) as fh:
        for l in fh:
            if l[0] == ">":
                if k is not None:
                    o[k] = "".join(b).upper()
                k = l[1:].split()[0]
                b = []
            else:
                b.append(l.strip())
    if k is not None:
        o[k] = "".join(b).upper()
    return o


def load_manifest_paths():
    m = {}
    with open(SRC + "/genome_map.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                m[f[0]] = f[1]
    return m


def event_delta():
    out = {}
    p = W + "/audits/delta_nonzero_classified.tsv"
    if os.path.exists(p):
        with open(p) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(cols)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out[f[ix["db_id"]]] = int(f[ix["delta"]])
    return out


def rows(key_col, shard, nshard):
    src = V2 + "/occurrences_flagged.tsv"
    by = collections.defaultdict(list)
    with open(src) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            k = f[ix[key_col]]
            if zlib.crc32(k.encode()) % nshard != shard % nshard:
                continue
            by[k].append(f)
    return ix, by


def seq_obj(seq_id, origin, s, segments, extra=None):
    o = {"sequence_id": seq_id, "origin": origin, "length_bp": len(s),
         "seq_sha256": sha(s), "segments": segments}
    if extra:
        o.update(extra)
    return o


def seg(acc, rec, a, b):
    return {"source_id": sid(acc, rec), "start0": a, "end0": b, "strand": "+"}


# --------------------------------------------------------------- pass carrier
def pass_carrier(a):
    man = load_manifest_paths()
    ix, by = rows("carrier_accession", a.shard, a.nshard)
    jl = open("%s/carrier_%03d.jsonl" % (a.out, a.shard), "w")
    fo = open("%s/carrier_%03d.fna" % (a.out, a.shard), "w")
    n = collections.Counter()
    accs = sorted(by)
    for gi, acc in enumerate(accs, 1):
        p = man.get(acc)
        ctgs = None
        if p:
            try:
                ctgs = fa(p)
            except (IOError, OSError):
                ctgs = None
        for f in by[acc]:
            occ = f[ix["occurrence_id"]]
            ctg = f[ix["carrier_contig"]]
            s, e = int(f[ix["insert_start"]]), int(f[ix["insert_end"]])
            ul, dl = int(f[ix["carrier_up_len"]]), int(f[ix["carrier_dn_len"]])
            rec = {"occurrence_id": occ, "sequences": {},
                   "missing_sequence_reasons": {}, "carrier_ok": False}
            S = ctgs.get(ctg) if ctgs else None
            if S is None:
                why = ("genome_not_in_manifest" if not p else
                       "contig_not_found_in_assembly")
                for r in ("insertion", "upstream", "downstream",
                          "reconstructed_empty"):
                    rec["sequences"][r] = None
                    rec["missing_sequence_reasons"][r] = why
                n[why] += 1
                jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
                continue
            L = len(S)
            if not (0 <= s < e <= L):
                for r in ("insertion", "upstream", "downstream",
                          "reconstructed_empty"):
                    rec["sequences"][r] = None
                    rec["missing_sequence_reasons"][r] = \
                        "insert_interval_outside_contig"
                n["insert_interval_outside_contig"] += 1
                jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
                continue
            us, ue = max(0, s - ul), s
            ds, de = e, min(L, e + dl)
            ins, up, dn = S[s:e], S[us:ue], S[ds:de]
            rec["carrier_ok"] = True
            rec["carrier_contig_len"] = L
            rec["intervals"] = {"up_start0": us, "up_end0": ue,
                                "down_start0": ds, "down_end0": de}
            rec["sequences"]["insertion"] = seq_obj(
                occ + ".insertion", "observed", ins, [seg(acc, ctg, s, e)])
            rec["sequences"]["upstream"] = seq_obj(
                occ + ".upstream", "observed", up, [seg(acc, ctg, us, ue)])
            rec["sequences"]["downstream"] = seq_obj(
                occ + ".downstream", "observed", dn, [seg(acc, ctg, ds, de)])
            rec["sequences"]["reconstructed_empty"] = seq_obj(
                occ + ".reconstructed_empty", "inferred", up + dn,
                [seg(acc, ctg, us, ue), seg(acc, ctg, ds, de)],
                {"native_locus": None, "construction_method": RECIPE,
                 "junction_offset0": len(up)})
            for r, sq in (("insertion", ins), ("upstream", up),
                          ("downstream", dn),
                          ("reconstructed_empty", up + dn)):
                fo.write(">%s.%s\n%s\n" % (occ, r, sq))
            n["ok"] += 1
            jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
        if gi % 200 == 0:
            sys.stderr.write("  %d/%d carriers\n" % (gi, len(accs)))
    jl.close()
    fo.close()
    print("carrier shard %d: assemblies %d" % (a.shard, len(accs)))
    for k, v in sorted(n.items()):
        print("   %-38s %7d" % (k, v))


# ------------------------------------------------------------- pass reference
def pass_reference(a):
    man = load_manifest_paths()
    ix, by = rows("empty_accession", a.shard, a.nshard)
    jl = open("%s/reference_%03d.jsonl" % (a.out, a.shard), "w")
    fo = open("%s/reference_%03d.fna" % (a.out, a.shard), "w")
    n = collections.Counter()
    accs = sorted(by)
    for gi, acc in enumerate(accs, 1):
        p = man.get(acc)
        ctgs = None
        if p:
            try:
                ctgs = fa(p)
            except (IOError, OSError):
                ctgs = None
        for f in by[acc]:
            occ = f[ix["occurrence_id"]]
            ctg = f[ix["empty_contig"]]
            j = int(f[ix["empty_junction_L"]])
            ul, dl = int(f[ix["carrier_up_len"]]), int(f[ix["carrier_dn_len"]])
            rec = {"occurrence_id": occ, "reference_empty": None,
                   "reason": None}
            S = ctgs.get(ctg) if ctgs else None
            if S is None:
                rec["reason"] = ("genome_not_in_manifest" if not p else
                                 "contig_not_found_in_assembly")
                n[rec["reason"]] += 1
                jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
                continue
            L = len(S)
            if not (0 <= j <= L):
                rec["reason"] = "junction_outside_contig"
                n["junction_outside_contig"] += 1
                jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
                continue
            ws, we = max(0, j - ul), min(L, j + dl)
            sq = S[ws:we]
            if not sq:
                rec["reason"] = "empty_reference_window"
                n["empty_reference_window"] += 1
                jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
                continue
            rec["reference_contig_len"] = L
            rec["window"] = {"start0": ws, "end0": we, "junction0": j}
            rec["reference_empty"] = seq_obj(
                occ + ".reference_empty", "observed", sq,
                [seg(acc, ctg, ws, we)], {"junction_offset0": j - ws})
            fo.write(">%s.reference_empty\n%s\n" % (occ, sq))
            n["ok"] += 1
            jl.write(json.dumps(rec, separators=(",", ":")) + "\n")
        if gi % 200 == 0:
            sys.stderr.write("  %d/%d references\n" % (gi, len(accs)))
    jl.close()
    fo.close()
    print("reference shard %d: assemblies %d" % (a.shard, len(accs)))
    for k, v in sorted(n.items()):
        print("   %-38s %7d" % (k, v))


# ----------------------------------------------------------------- pass merge
VIEW = ["occurrence_id", "event_id", "analysis_set", "official",
        "carrier_assembly_accession_version", "carrier_sequence_accession_version",
        "carrier_fasta_record_id", "insert_start0", "insert_end0",
        "insertion_sequence_id", "insertion_seq",
        "up_start0", "up_end0", "upstream_sequence_id", "upstream_seq",
        "down_start0", "down_end0", "downstream_sequence_id", "downstream_seq",
        "reference_assembly_accession_version",
        "reference_sequence_accession_version", "reference_fasta_record_id",
        "reference_window_start0", "reference_window_end0",
        "reference_junction0", "reference_to_carrier_strand",
        "reference_empty_sequence_id", "reference_empty_seq",
        "reconstructed_empty_sequence_id", "reconstructed_empty_seq",
        "reconstructed_junction_offset0", "reconstructed_source_segments_json",
        "breakpoint_status", "min_joint_shift_bp", "max_joint_shift_bp",
        "chosen_joint_shift_bp", "delta_band", "auditA_verdict",
        "exclusion_reasons_json", "source_validation_status",
        "recipe_validation_status"]
NOSEQ = {"insertion_seq", "upstream_seq", "downstream_seq",
         "reference_empty_seq", "reconstructed_empty_seq"}


def pass_merge(a):
    delta = event_delta()
    ref = {}
    for p in sorted(os.listdir(a.out)):
        if p.startswith("reference_") and p.endswith(".jsonl"):
            with open(a.out + "/" + p) as fh:
                for line in fh:
                    d = json.loads(line)
                    ref[d["occurrence_id"]] = d
    print("reference fragments loaded: %d" % len(ref))

    seqs = {}
    with open(V2 + "/occurrences_flagged.tsv") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            seqs[f[ix["occurrence_id"]]] = f
    print("input rows: %d" % len(seqs))

    body = {}
    for p in sorted(os.listdir(a.out)):
        if p.startswith("carrier_") and p.endswith(".jsonl"):
            with open(a.out + "/" + p) as fh:
                for line in fh:
                    d = json.loads(line)
                    body[d["occurrence_id"]] = d
    print("carrier fragments loaded: %d" % len(body))

    txt = {}
    if a.inline_seq:
        for p in sorted(os.listdir(a.out)):
            if p.endswith(".fna") and (p.startswith("carrier_") or
                                       p.startswith("reference_")):
                k = None
                with open(a.out + "/" + p) as fh:
                    for line in fh:
                        if line[0] == ">":
                            k = line[1:].strip()
                        else:
                            txt[k] = line.strip()

    oj = open(OUT + "/occurrences.jsonl", "w")
    vw = open(OUT + "/occurrences_view.tsv", "w")
    br = open(OUT + "/occurrences_browse.tsv", "w")
    re_ev = open(OUT + "/reference_evidence.jsonl", "w")
    eq_ev = open(OUT + "/equivalence_evidence.jsonl", "w")
    vw.write("\t".join(VIEW) + "\n")
    bcols = [c for c in VIEW if c not in NOSEQ]
    br.write("\t".join(bcols) + "\n")

    n = collections.Counter()
    nseq = collections.Counter()
    for occ in sorted(seqs):
        f = seqs[occ]
        c = body.get(occ)
        r = ref.get(occ)
        if c is None:
            n["MISSING_CARRIER_FRAGMENT"] += 1
            continue
        n["records"] += 1
        acc, ctg = f[ix["carrier_accession"]], f[ix["carrier_contig"]]
        racc, rctg = f[ix["empty_accession"]], f[ix["empty_contig"]]
        ev_ref = "REFEV|" + occ
        ev_eq = "EQEV|" + occ
        bstat = f[ix["breakpoint_status"]]
        el, er = int(f[ix["carrier_equiv_left"]]), int(f[ix["carrier_equiv_right"]])
        if bstat == "unique":
            mn, mx = 0, 0
        elif bstat in ("equivalent_range", "beyond_equivalence_range",
                       "discordant_sides"):
            mn, mx = el, er
        else:
            mn, mx = None, None
        v2shift = f[ix["chosen_shift"]]
        sq = dict(c["sequences"])
        miss = dict(c["missing_sequence_reasons"])
        if r is None:
            sq["reference_empty"] = None
            miss["reference_empty"] = "reference_fragment_absent"
        elif r.get("reference_empty") is None:
            sq["reference_empty"] = None
            miss["reference_empty"] = r.get("reason") or "unknown"
        else:
            sq["reference_empty"] = r["reference_empty"]
        for role in ROLES:
            nseq["present" if sq.get(role) else "null"] += 1

        rec = {
            "schema_version": SCHEMA,
            "occurrence_id": occ,
            "event_id": f[ix["db_id"]],
            "record_revision": 1,
            "carrier": {
                "source_id": sid(acc, ctg),
                "original_insert_start0": int(f[ix["insert_start"]]),
                "original_insert_end0": int(f[ix["insert_end"]]),
                "insert_start0": int(f[ix["insert_start"]]),
                "insert_end0": int(f[ix["insert_end"]]),
                "canonical_insert_strand_in_carrier":
                    f[ix["insert_strand"]] if f[ix["insert_strand"]] in "+-"
                    else "unknown",
                "interval_definition": INTERVAL_DEF,
            },
            "reference": {
                "source_id": sid(racc, rctg),
                "junction0": int(f[ix["empty_junction_L"]]),
                "reference_to_carrier_strand":
                    {"same": "+", "flip": "-"}.get(f[ix["orientation"]],
                                                   "unknown"),
                "selection_evidence_id": ev_ref,
            },
            "breakpoint": {
                "status": bstat if bstat else "unresolved",
                "min_joint_shift_bp": mn,
                "max_joint_shift_bp": mx,
                # ALWAYS 0: the exported bases were cut at the original
                # interval. See the module docstring.
                "chosen_joint_shift_bp": 0,
                "v2_reported_equivalent_shift_bp":
                    (int(v2shift) if v2shift not in (".", "") else None),
                "equivalence_evidence_id": ev_eq,
            },
            "sequences": sq,
            "quality": {
                "official": f[ix["official_fasta"]] == "yes",
                "analysis_set": f[ix["analysis_set"]],
                "localization_status": f[ix["localization_status"]],
                "reconstructed_sequence_status":
                    f[ix["reconstructed_sequence_status"]],
                "candidate_search_truncated":
                    f[ix["candidate_search_truncated"]] == "yes",
                "core_recoverable":
                    {"yes": True, "no": False}.get(f[ix["core_recoverable"]],
                                                   "unknown"),
                "delta_band": f[ix["delta_band"]],
                "event_delta_bp": delta.get(f[ix["db_id"]], 0),
                "auditA_verdict": f[ix["auditA_verdict"]],
                "exclusion_reasons":
                    [] if f[ix["exclusion_reason"]] == "." else
                    f[ix["exclusion_reason"]].split(";"),
            },
            "missing_sequence_reasons": miss,
        }
        oj.write(json.dumps(rec, separators=(",", ":")) + "\n")
        re_ev.write(json.dumps({
            "evidence_id": ev_ref, "occurrence_id": occ,
            "n_candidates_examined": f[ix["n_candidates_examined"]],
            "n_candidates_trusted": f[ix["n_candidates_trusted"]],
            "candidate_search_truncated": f[ix["candidate_search_truncated"]],
            "orientation": f[ix["orientation"]], "anchor_bp": f[ix["anchor_bp"]],
            "up_identity": f[ix["up_identity"]], "up_coverage": f[ix["up_coverage"]],
            "dn_identity": f[ix["dn_identity"]], "dn_coverage": f[ix["dn_coverage"]],
            "margin": f[ix["margin"]],
            "localization_status": f[ix["localization_status"]],
        }, separators=(",", ":")) + "\n")
        eq_ev.write(json.dumps({
            "evidence_id": ev_eq, "occurrence_id": occ,
            "carrier_equiv_left": el, "carrier_equiv_right": er,
            "up_shift": f[ix["up_shift"]], "dn_shift": f[ix["dn_shift"]],
            "shift_coherent": f[ix["shift_coherent"]],
            "sequence_equivalence_verified":
                f[ix["sequence_equivalence_verified"]],
            "breakpoint_status": bstat,
        }, separators=(",", ":")) + "\n")

        iv = c.get("intervals") or {}
        w = (r or {}).get("window") or {}
        row = {
            "occurrence_id": occ, "event_id": f[ix["db_id"]],
            "analysis_set": f[ix["analysis_set"]],
            "official": f[ix["official_fasta"]],
            "carrier_assembly_accession_version": acc,
            "carrier_sequence_accession_version": "",
            "carrier_fasta_record_id": ctg,
            "insert_start0": f[ix["insert_start"]],
            "insert_end0": f[ix["insert_end"]],
            "up_start0": iv.get("up_start0", ""), "up_end0": iv.get("up_end0", ""),
            "down_start0": iv.get("down_start0", ""),
            "down_end0": iv.get("down_end0", ""),
            "reference_assembly_accession_version": racc,
            "reference_sequence_accession_version": "",
            "reference_fasta_record_id": rctg,
            "reference_window_start0": w.get("start0", ""),
            "reference_window_end0": w.get("end0", ""),
            "reference_junction0": w.get("junction0", f[ix["empty_junction_L"]]),
            "reference_to_carrier_strand":
                rec["reference"]["reference_to_carrier_strand"],
            "reconstructed_junction_offset0":
                (sq["reconstructed_empty"] or {}).get("junction_offset0", ""),
            "reconstructed_source_segments_json":
                json.dumps((sq["reconstructed_empty"] or {}).get("segments", []),
                           separators=(",", ":")),
            "breakpoint_status": bstat,
            "min_joint_shift_bp": "" if mn is None else mn,
            "max_joint_shift_bp": "" if mx is None else mx,
            "chosen_joint_shift_bp": 0,
            "delta_band": f[ix["delta_band"]],
            "auditA_verdict": f[ix["auditA_verdict"]],
            "exclusion_reasons_json":
                json.dumps(rec["quality"]["exclusion_reasons"],
                           separators=(",", ":")),
            "source_validation_status": "NOT_RUN",
            "recipe_validation_status": "NOT_RUN",
        }
        for role, pre in (("insertion", "insertion"), ("upstream", "upstream"),
                          ("downstream", "downstream"),
                          ("reference_empty", "reference_empty"),
                          ("reconstructed_empty", "reconstructed_empty")):
            o = sq.get(role)
            row[pre + "_sequence_id"] = o["sequence_id"] if o else ""
            row[pre + "_seq"] = (txt.get(o["sequence_id"], "") if
                                 (o and a.inline_seq) else "")
        vw.write("\t".join(str(row.get(k, "")) for k in VIEW) + "\n")
        br.write("\t".join(str(row.get(k, "")) for k in bcols) + "\n")
    for h in (oj, vw, br, re_ev, eq_ev):
        h.close()
    print("merged records: %d" % n["records"])
    if n["MISSING_CARRIER_FRAGMENT"]:
        print("   MISSING_CARRIER_FRAGMENT %d" % n["MISSING_CARRIER_FRAGMENT"])
    print("   sequence objects present %d  null %d"
          % (nseq["present"], nseq["null"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="phase", required=True,
                    choices=("carrier", "reference", "merge"))
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--inline-seq", action="store_true",
                    help="fill the *_seq columns of occurrences_view.tsv")
    a = ap.parse_args()
    {"carrier": pass_carrier, "reference": pass_reference,
     "merge": pass_merge}[a.phase](a)


if __name__ == "__main__":
    main()
