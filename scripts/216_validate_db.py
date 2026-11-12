#!/usr/bin/env python3
"""Stage 3: acceptance for the v1 record-format export (spec section 8).

Passes:
  --pass sources   B. frozen sources: file byte hash, contig hash, contig
                   length, duplicate record ids. Each file hashed once.
  --pass records   A/C/D/E per occurrence against the ORIGINAL FNA.
  --pass summarise F. merge shard results, write validation.tsv and
                   validation_summary.json, exit non-zero on any FAIL.

WHAT THIS PROVES AND WHAT IT DOES NOT. Every PASS here says one thing: the
stored bases are the bases at the declared coordinates of the declared
assembly. It says nothing about whether the empty/insertion pairing is
biologically right. Those stay separate columns and separate totals, and the
source-validation rate is never reported as an insertion accuracy.

A missing or unsupported check is NOT_VERIFIABLE, never PASS. Records whose
sequence objects are null still count toward coverage.
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
RECIPE = "concat_carrier_flanks_v2"
OBS = ("insertion", "upstream", "downstream", "reference_empty")
ROLES = OBS + ("reconstructed_empty",)
REPLAY_CAP = 200      # legal shifts replayed per record, matching 206's cap
REPLAY_FLANK = 300    # fixed outer window used for every k of one record


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def fa_strict(p):
    out, k, b, dup = {}, None, [], []
    with open(p) as fh:
        for l in fh:
            if l[0] == ">":
                if k is not None:
                    if k in out:
                        dup.append(k)
                    out[k] = "".join(b).upper()
                k = l[1:].split()[0]
                b = []
            else:
                b.append(l.strip())
    if k is not None:
        if k in out:
            dup.append(k)
        out[k] = "".join(b).upper()
    return out, dup


def file_sha256(p, buf=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            x = fh.read(buf)
            if not x:
                break
            h.update(x)
    return h.hexdigest()


def load_sources():
    m = {}
    with open(OUT + "/source_manifest.jsonl") as fh:
        for line in fh:
            d = json.loads(line)
            m[d["source_id"]] = d
    return m


# --------------------------------------------------------------- pass sources
def pass_sources(a):
    srcs = load_sources()
    byfile = collections.defaultdict(list)
    for d in srcs.values():
        byfile[d["fna_path"]].append(d)
    paths = sorted(p for p in byfile
                   if zlib.crc32(p.encode()) % a.nshard == a.shard % a.nshard)
    out = open("%s/src_%03d.tsv" % (a.out, a.shard), "w")
    n = collections.Counter()
    for i, rel in enumerate(paths, 1):
        full = os.path.join(W, rel)
        if not os.path.exists(full):
            for d in byfile[rel]:
                out.write("%s\tfile\tFAIL\tfile_missing\n" % d["source_id"])
                n["FAIL"] += 1
            continue
        fh256 = file_sha256(full)
        ctgs, dup = fa_strict(full)
        for d in byfile[rel]:
            sidv, rec = d["source_id"], d["fasta_record_id"]
            errs = []
            if fh256 != d["fna_file_sha256"]:
                errs.append("fna_file_sha256_mismatch")
            if rec in dup:
                errs.append("duplicate_fasta_record_id")
            s = ctgs.get(rec)
            if s is None:
                errs.append("fasta_record_id_not_found")
            else:
                if len(s) != d["contig_length_bp"]:
                    errs.append("contig_length_mismatch")
                if sha(s) != d["contig_seq_sha256"]:
                    errs.append("contig_seq_sha256_mismatch")
            st = "FAIL" if errs else "PASS"
            n[st] += 1
            out.write("%s\tsource\t%s\t%s\n" % (sidv, st, ",".join(errs) or "."))
        if i % 200 == 0:
            sys.stderr.write("  %d/%d files\n" % (i, len(paths)))
    out.close()
    print("sources shard %d: files %d  PASS %d  FAIL %d"
          % (a.shard, len(paths), n["PASS"], n["FAIL"]))


# --------------------------------------------------------------- pass records
def replay(S, s, e, mn, mx, flank=REPLAY_FLANK, cap=REPLAY_CAP):
    """Replay every legal k in [mn,mx] in ONE fixed outer window.

    The window must not move with k, or the comparison is between two
    different excerpts rather than two different cuts of the same excerpt.
    Returns (all_equal, n_tested, range_complete).
    """
    lo, hi = s - flank, e + flank
    if lo < 0 or hi > len(S):
        return None, 0, False
    base = S[lo:s] + S[e:hi]
    ks = [k for k in range(mn, mx + 1)]
    complete = True
    if len(ks) > cap:
        ks, complete = ks[:cap], False
    ok, n = True, 0
    for k in ks:
        if s + k < lo or e + k > hi:
            complete = False
            continue
        if S[lo:s + k] + S[e + k:hi] != base:
            ok = False
        n += 1
    return ok, n, complete


def pass_records(a):
    srcs = load_sources()
    man = {}
    with open(SRC + "/genome_map.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                man[f[0]] = f[1]

    # Shard by the source this pass actually opens. Sharding the reference
    # checks by CARRIER would scatter the reference assemblies across every
    # shard and reload most of them 16 times.
    refs_only = a.phase == "refrecords"

    def key(d):
        if refs_only:
            o = d["sequences"].get("reference_empty")
            return o["segments"][0]["source_id"] if o else None
        return d["carrier"]["source_id"]

    recs = {}
    with open(OUT + "/occurrences.jsonl") as fh:
        for line in fh:
            d = json.loads(line)
            k = key(d)
            if k is None:
                continue
            if zlib.crc32(k.encode()) % a.nshard == a.shard % a.nshard:
                recs[d["occurrence_id"]] = d
    want = set()
    for d in recs.values():
        for r in (("reference_empty",) if refs_only else ROLES):
            o = d["sequences"].get(r)
            if o:
                want.add(o["sequence_id"])
    txt, k = {}, None
    with open(OUT + "/sequences.fna") as fh:
        for line in fh:
            if line[0] == ">":
                k = line[1:].strip()
                k = k if k in want else None
            elif k:
                txt[k] = line.strip()
    sys.stderr.write("shard %d: %d records, %d sequences loaded\n"
                     % (a.shard, len(recs), len(txt)))

    bycar = collections.defaultdict(list)
    byref = collections.defaultdict(list)
    for d in recs.values():
        if refs_only:
            o = d["sequences"].get("reference_empty")
            byref[o["segments"][0]["source_id"]].append(d)
        else:
            bycar[d["carrier"]["source_id"]].append(d)

    tag = "ref" if refs_only else "car"
    out = open("%s/val_%s_%03d.tsv" % (a.out, tag, a.shard), "w")
    summ = open("%s/rec_%s_%03d.tsv" % (a.out, tag, a.shard), "w")

    def emit(occ, role, kind, status, err, sidv, st, en, elen, alen, eh, ah):
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                  % (occ, role, kind, status, err, sidv, st, en, elen, alen,
                     eh, ah))

    def slice_check(d, role, S, clen):
        o = d["sequences"].get(role)
        occ = d["occurrence_id"]
        if o is None:
            emit(occ, role, "observed_slice", "NOT_VERIFIABLE",
                 d["missing_sequence_reasons"].get(role, "null_sequence"),
                 ".", ".", ".", ".", ".", ".", ".")
            return None
        sg = o["segments"][0]
        st, en = sg["start0"], sg["end0"]
        errs = []
        if not (0 <= st < en <= clen):
            errs.append("interval_outside_contig")
            got = ""
        else:
            got = S[st:en]
        exp = txt.get(o["sequence_id"])
        if exp is None:
            errs.append("sequence_id_absent_from_fasta")
            exp = ""
        if got != exp:
            errs.append("sequence_differs_from_declared_coordinates")
        if len(exp) != o["length_bp"]:
            errs.append("length_mismatch")
        if exp and sha(exp) != o["seq_sha256"]:
            errs.append("sha256_mismatch")
        st_ = "FAIL" if errs else "PASS"
        emit(occ, role, "observed_slice", st_, ",".join(errs) or ".",
             sg["source_id"], st, en, o["length_bp"], len(exp),
             o["seq_sha256"], sha(exp) if exp else ".")
        return None if errs else got

    n = collections.Counter()
    for sidv, ds in bycar.items():
        meta = srcs.get(sidv)
        S, clen = None, 0
        if meta:
            p = os.path.join(W, meta["fna_path"])
            try:
                ctgs, _ = fa_strict(p)
                S = ctgs.get(meta["fasta_record_id"])
                clen = len(S) if S is not None else 0
            except (IOError, OSError):
                S = None
        for d in ds:
            occ = d["occurrence_id"]
            res = {"observed_slice_pass": 0, "coordinate_geometry_pass": 0,
                   "inferred_recipe_pass": 0, "breakpoint_replay_pass": 0,
                   "not_verifiable": 0, "fail": 0}
            if S is None:
                # Only the carrier-derived roles; reference_empty belongs to
                # the reference pass and would otherwise be emitted twice.
                for r in ("insertion", "upstream", "downstream",
                          "reconstructed_empty"):
                    emit(occ, r, "observed_slice", "NOT_VERIFIABLE",
                         "carrier_source_unreadable", sidv, ".", ".", ".",
                         ".", ".", ".")
                    res["not_verifiable"] += 1
                summ.write("%s\t%s\n" % (occ, json.dumps(res)))
                n["carrier_unreadable"] += 1
                continue
            ins = slice_check(d, "insertion", S, clen)
            up = slice_check(d, "upstream", S, clen)
            dn = slice_check(d, "downstream", S, clen)
            for v in (ins, up, dn):
                res["observed_slice_pass" if v is not None else "fail"] += 1
            # ---- D. adjacency + recipe, re-read from the FNA ----
            o = d["sequences"].get("reconstructed_empty")
            errs = []
            if o is None:
                emit(occ, "reconstructed_empty", "inferred_recipe",
                     "NOT_VERIFIABLE",
                     d["missing_sequence_reasons"].get(
                         "reconstructed_empty", "null_sequence"),
                     ".", ".", ".", ".", ".", ".", ".")
                res["not_verifiable"] += 1
            else:
                segs = o["segments"]
                c = d["carrier"]
                if len(segs) != 2:
                    errs.append("recipe_segment_count")
                else:
                    a1, a2 = segs
                    if a1["source_id"] != sidv or a2["source_id"] != sidv:
                        errs.append("recipe_segments_not_same_carrier")
                    if a1["end0"] != c["insert_start0"]:
                        errs.append("upstream_not_adjacent_to_insert_start")
                    if a2["start0"] != c["insert_end0"]:
                        errs.append("downstream_not_adjacent_to_insert_end")
                    reb = S[a1["start0"]:a1["end0"]] + S[a2["start0"]:a2["end0"]]
                    if reb != txt.get(o["sequence_id"], ""):
                        errs.append("recipe_replay_differs")
                    if len(reb) != o["length_bp"]:
                        errs.append("length_mismatch")
                    if sha(reb) != o["seq_sha256"]:
                        errs.append("sha256_mismatch")
                    if o.get("junction_offset0") != a1["end0"] - a1["start0"]:
                        errs.append("junction_offset_mismatch")
                    whole = S[a1["start0"]:a2["end0"]]
                    if whole != (S[a1["start0"]:a1["end0"]] +
                                 S[c["insert_start0"]:c["insert_end0"]] +
                                 S[a2["start0"]:a2["end0"]]):
                        errs.append("up_insert_down_not_contiguous")
                if o.get("native_locus", "MISSING") is not None:
                    errs.append("inferred_sequence_claims_native_locus")
                if o.get("construction_method") != RECIPE:
                    errs.append("construction_method_mismatch")
                if o.get("origin") != "inferred":
                    errs.append("origin_not_inferred")
                st_ = "FAIL" if errs else "PASS"
                emit(occ, "reconstructed_empty", "inferred_recipe", st_,
                     ",".join(errs) or ".", sidv, ".", ".", o["length_bp"],
                     ".", o["seq_sha256"], ".")
                res["inferred_recipe_pass" if not errs else "fail"] += 1
            # ---- geometry ----
            c = d["carrier"]
            g = []
            if c["insert_start0"] != c["original_insert_start0"] + \
                    d["breakpoint"]["chosen_joint_shift_bp"]:
                g.append("start_not_original_plus_shift")
            if c["insert_end0"] != c["original_insert_end0"] + \
                    d["breakpoint"]["chosen_joint_shift_bp"]:
                g.append("end_not_original_plus_shift")
            if not (0 <= c["insert_start0"] < c["insert_end0"] <= clen):
                g.append("insert_interval_outside_contig")
            emit(occ, "carrier", "coordinate_geometry",
                 "FAIL" if g else "PASS", ",".join(g) or ".", sidv,
                 c["insert_start0"], c["insert_end0"], ".", ".", ".", ".")
            res["coordinate_geometry_pass" if not g else "fail"] += 1
            # ---- E. breakpoint replay ----
            b = d["breakpoint"]
            mn, mx = b["min_joint_shift_bp"], b["max_joint_shift_bp"]
            if mn is None or mx is None:
                emit(occ, "carrier", "breakpoint_replay", "NOT_VERIFIABLE",
                     "no_declared_range", sidv, ".", ".", ".", ".", ".", ".")
                res["not_verifiable"] += 1
            else:
                ok, ntest, comp = replay(S, c["insert_start0"],
                                         c["insert_end0"], mn, mx,
                                         flank=a.replay_flank)
                if ok is None:
                    emit(occ, "carrier", "breakpoint_replay", "NOT_VERIFIABLE",
                         "replay_window_clipped_by_contig", sidv, mn, mx,
                         ".", ".", ".", ".")
                    res["not_verifiable"] += 1
                else:
                    e2 = [] if ok else ["declared_range_not_sequence_equivalent"]
                    if not comp:
                        e2.append("range_complete=false")
                    emit(occ, "carrier", "breakpoint_replay",
                         "FAIL" if not ok else "PASS", ",".join(e2) or ".",
                         sidv, mn, mx, ntest, ".", ".", ".")
                    res["breakpoint_replay_pass" if ok else "fail"] += 1
            summ.write("%s\t%s\n" % (occ, json.dumps(res)))
            n["records"] += 1

    # ---- reference_empty slices, grouped by reference source ----
    for sidv, ds in byref.items():
        meta = srcs.get(sidv)
        S, clen = None, 0
        if meta:
            p = os.path.join(W, meta["fna_path"])
            try:
                ctgs, _ = fa_strict(p)
                S = ctgs.get(meta["fasta_record_id"])
                clen = len(S) if S is not None else 0
            except (IOError, OSError):
                S = None
        for d in ds:
            occ = d["occurrence_id"]
            if S is None:
                emit(occ, "reference_empty", "observed_slice",
                     "NOT_VERIFIABLE", "reference_source_unreadable", sidv,
                     ".", ".", ".", ".", ".", ".")
                continue
            got = slice_check(d, "reference_empty", S, clen)
            o = d["sequences"]["reference_empty"]
            sg = o["segments"][0]
            j = d["reference"]["junction0"]
            e3 = []
            if o.get("junction_offset0") != j - sg["start0"]:
                e3.append("reference_junction_offset_mismatch")
            if not (sg["start0"] <= j <= sg["end0"]):
                e3.append("junction_outside_window")
            emit(occ, "reference_empty", "reference_junction_geometry",
                 "FAIL" if e3 else "PASS", ",".join(e3) or ".", sidv,
                 sg["start0"], sg["end0"], j, o.get("junction_offset0"),
                 ".", ".")
    out.close()
    summ.close()
    print("records shard %d: %d records, %d carriers unreadable"
          % (a.shard, n["records"], n["carrier_unreadable"]))


# ------------------------------------------------------------- pass summarise
def pass_summarise(a):
    expected, official, conservative = set(), set(), set()
    with open(V2 + "/occurrences_flagged.tsv") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            o = f[ix["occurrence_id"]]
            expected.add(o)
            if f[ix["official_fasta"]] == "yes":
                official.add(o)
            if f[ix["analysis_set"]] == "conservative":
                conservative.add(o)

    got, dup = set(), 0
    refd = set()
    qoff, qcons = set(), set()
    with open(OUT + "/occurrences.jsonl") as fh:
        for line in fh:
            d = json.loads(line)
            o = d["occurrence_id"]
            if o in got:
                dup += 1
            got.add(o)
            if d["quality"]["official"]:
                qoff.add(o)
            if d["quality"]["analysis_set"] == "conservative":
                qcons.add(o)
            for r in ROLES:
                x = d["sequences"].get(r)
                if x:
                    refd.add(x["sequence_id"])

    fids, fdup = set(), 0
    with open(OUT + "/sequences.fna") as fh:
        for line in fh:
            if line[0] == ">":
                k = line[1:].strip()
                if k in fids:
                    fdup += 1
                fids.add(k)

    print("A. COVERAGE AND IDENTITY (set relations, not counts)")
    print("   expected input ids      %7d" % len(expected))
    print("   records emitted         %7d   duplicates %d" % (len(got), dup))
    print("   missing                 %7d" % len(expected - got))
    print("   unexpected              %7d" % len(got - expected))
    print("   official  expected %d  in records %d  symmetric diff %d"
          % (len(official), len(qoff), len(official ^ qoff)))
    print("   conservative expected %d  in records %d  symmetric diff %d"
          % (len(conservative), len(qcons), len(conservative ^ qcons)))
    print("   sequences.fna ids       %7d   duplicates %d" % (len(fids), fdup))
    print("   referenced but absent   %7d" % len(refd - fids))
    print("   present but unreferenced%7d" % len(fids - refd))

    st = collections.Counter()
    bykind = collections.Counter()
    errs = collections.Counter()
    failrows = []
    vt = open(OUT + "/validation.tsv", "w")
    vt.write("occurrence_id\trole\tvalidation_kind\tstatus\terror_code\t"
             "source_id\tstart0\tend0\texpected_len\tactual_len\t"
             "expected_sha256\tactual_sha256\n")
    for p in sorted(os.listdir(a.out)):
        if not (p.startswith("val_") and p.endswith(".tsv")):
            continue
        with open(a.out + "/" + p) as fh:
            for line in fh:
                vt.write(line)
                f = line.rstrip("\n").split("\t")
                st[f[3]] += 1
                bykind[(f[2], f[3])] += 1
                if f[3] == "FAIL":
                    errs[f[4]] += 1
                    if len(failrows) < 50:
                        failrows.append(f[:5])
    for p in sorted(os.listdir(a.out)):
        if p.startswith("src_") and p.endswith(".tsv"):
            with open(a.out + "/" + p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    vt.write("%s\t.\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\n"
                             % (f[0], f[1], f[2], f[3], f[0]))
                    st[f[2]] += 1
                    bykind[("source", f[2])] += 1
                    if f[2] == "FAIL":
                        errs[f[3]] += 1
    vt.close()

    print("\nC-F. PER-CHECK RESULTS")
    for (k, s), v in sorted(bykind.items()):
        print("   %-32s %-16s %8d" % (k, s, v))
    print("\n   TOTAL  PASS %d  FAIL %d  NOT_VERIFIABLE %d"
          % (st["PASS"], st["FAIL"], st["NOT_VERIFIABLE"]))
    if errs:
        print("\n   failure codes")
        for k, v in errs.most_common(20):
            print("      %-54s %7d" % (k, v))
    for f in failrows[:20]:
        print("      FAIL %s" % "\t".join(f))

    summary = {
        "scope": "full_corpus",
        "biological_accuracy_validated": False,
        "source_validation_is_not_insertion_accuracy": True,
        "expected_records": len(expected),
        "emitted_records": len(got),
        "missing_records": len(expected - got),
        "unexpected_records": len(got - expected),
        "duplicate_records": dup,
        "official_symmetric_difference": len(official ^ qoff),
        "conservative_symmetric_difference": len(conservative ^ qcons),
        "sequence_ids": len(fids),
        "duplicate_sequence_ids": fdup,
        "referenced_but_absent": len(refd - fids),
        "present_but_unreferenced": len(fids - refd),
        "pass": st["PASS"], "fail": st["FAIL"],
        "not_verifiable": st["NOT_VERIFIABLE"],
        "by_kind": {"%s|%s" % k: v for k, v in sorted(bykind.items())},
        "failure_codes": dict(errs),
    }
    with open(OUT + "/validation_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    bad = (st["FAIL"] or dup or fdup or (expected - got) or (got - expected)
           or (refd - fids) or (official ^ qoff) or (conservative ^ qcons))
    print("\nACCEPTANCE: %s" % ("FAILED" if bad else "PASSED"))
    if bad:
        sys.exit(1)


def main():
    global W, OUT, V2, SRC
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="phase", required=True,
                    choices=("sources", "records", "refrecords", "summarise"))
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--w", default=W, help="corpus root holding fna_path")
    ap.add_argument("--db", default=OUT, help="export package directory")
    ap.add_argument("--v2", default=V2, help="dir holding occurrences_flagged.tsv")
    ap.add_argument("--src", default=SRC, help="dir holding genome_map.tsv")
    ap.add_argument("--replay-flank", type=int, default=REPLAY_FLANK,
                    help="fixed outer window used for every k of one record")
    a = ap.parse_args()
    W, OUT, V2, SRC = a.w, a.db, a.v2, a.src
    {"sources": pass_sources, "records": pass_records,
     "refrecords": pass_records,
     "summarise": pass_summarise}[a.phase](a)


if __name__ == "__main__":
    main()
