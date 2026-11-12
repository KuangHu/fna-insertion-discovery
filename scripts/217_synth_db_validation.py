#!/usr/bin/env python3
"""Negative controls for the export validator (spec section 10's test list).

A validator that only ever sees correct input proves nothing. Each case below
builds a deliberately wrong package and asserts the validator REJECTS it, and
the two "allowed difference" cases assert it does NOT reject them:

  reverse canonical insert      must PASS  -- stored bases stay carrier-forward
  reverse reference alignment   must PASS  -- reference stays its own forward
  real SNP vs the reference     must PASS  -- cross-assembly identity is not required
  window offset wrong           must FAIL
  wrong assembly, same contig   must FAIL  -- the token is not an identifier
  coordinates off by one        must FAIL
  duplicate FASTA record id     must FAIL
  join order reversed           must FAIL
  invented native locus         must FAIL
  non-equivalent claimed range  must FAIL
  missing sequence              must be NOT_VERIFIABLE, never PASS
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
VAL = os.path.join(HERE, "216_validate_db.py")
RECIPE = "concat_carrier_flanks_v2"


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def fsha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


UP = "ACGTTGCAAGGCTTACCGATAGCTTAGGCCATTGACCTAGGATTCCAGTA"
INS = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG"
DN = "GGCATTACGTAGCTAGCTTAGCCATGGATCCGGTACGATTCAGGCTAACT"
# one SNP in the reference copy of the right flank
DN_SNP = DN[:10] + ("A" if DN[10] != "A" else "C") + DN[11:]


def build(root):
    os.makedirs(root + "/fna", exist_ok=True)
    os.makedirs(root + "/db", exist_ok=True)
    os.makedirs(root + "/v2", exist_ok=True)
    os.makedirs(root + "/src", exist_ok=True)
    os.makedirs(root + "/work", exist_ok=True)
    car = UP + INS + DN
    ref = UP + DN_SNP
    with open(root + "/fna/carrierA.fna", "w") as fh:
        fh.write(">contig1\n%s\n" % car)
    with open(root + "/fna/refB.fna", "w") as fh:
        fh.write(">contig1\n%s\n" % ref)          # SAME token, other assembly
    with open(root + "/fna/dupC.fna", "w") as fh:
        fh.write(">contig1\n%s\n>contig1\n%s\n" % (car, car))
    with open(root + "/src/genome_map.tsv", "w") as fh:
        fh.write("accession\tpath\n")
        for a, f in (("ASM_A", "carrierA.fna"), ("ASM_B", "refB.fna"),
                     ("ASM_C", "dupC.fna")):
            fh.write("%s\t%s/fna/%s\n" % (a, root, f))
    srcs = []
    for a, f, token, seq in (("ASM_A", "carrierA.fna", "contig1", car),
                             ("ASM_B", "refB.fna", "contig1", ref),
                             ("ASM_C", "dupC.fna", "contig1", car)):
        srcs.append({
            "source_id": "SRC|%s|%s" % (a, token),
            "assembly_accession_version": a, "assembly_local_id": None,
            "sequence_accession_version": None, "fna_path": "fna/" + f,
            "fasta_record_id": token, "contig_length_bp": len(seq),
            "contig_seq_sha256": sha(seq),
            "fna_file_sha256": fsha("%s/fna/%s" % (root, f)),
            "topology": "linear"})
    with open(root + "/db/source_manifest.jsonl", "w") as fh:
        for d in srcs:
            fh.write(json.dumps(d) + "\n")
    return car, ref


S0, E0 = len(UP), len(UP) + len(INS)
US, UE, DS = 0, len(UP), len(UP) + len(INS)


def obj(sid, origin, s, segs, extra=None):
    o = {"sequence_id": sid, "origin": origin, "length_bp": len(s),
         "seq_sha256": sha(s), "segments": segs}
    if extra:
        o.update(extra)
    return o


def seg(src, a, b):
    return {"source_id": src, "start0": a, "end0": b, "strand": "+"}


def record(occ, car, ref, **kw):
    A = "SRC|ASM_A|contig1"
    B = "SRC|ASM_B|contig1"
    carrier_src = kw.get("carrier_src", A)
    s = kw.get("ins_start", S0)
    e = kw.get("ins_end", E0)
    de = len(car)
    ins, up, dn = car[S0:E0], car[US:UE], car[DS:de]
    rsegs = [seg(carrier_src, US, UE), seg(carrier_src, DS, de)]
    if kw.get("reverse_join"):
        rsegs = rsegs[::-1]
    refw_s, refw_e, j = 0, len(ref), len(UP)
    joff = j - refw_s + (1 if kw.get("bad_junction_offset") else 0)
    seqs = {
        "insertion": obj(occ + ".insertion", "observed", ins,
                         [seg(carrier_src, s, e)]),
        "upstream": obj(occ + ".upstream", "observed", up,
                        [seg(carrier_src, US, UE)]),
        "downstream": obj(occ + ".downstream", "observed", dn,
                          [seg(carrier_src, DS, de)]),
        "reference_empty": obj(occ + ".reference_empty", "observed",
                               ref[refw_s:refw_e], [seg(B, refw_s, refw_e)],
                               {"junction_offset0": joff}),
        "reconstructed_empty": obj(
            occ + ".reconstructed_empty", "inferred", up + dn, rsegs,
            {"native_locus": kw.get("native_locus"),
             "construction_method": RECIPE, "junction_offset0": len(up)}),
    }
    miss = {}
    if kw.get("drop_insertion"):
        seqs["insertion"] = None
        miss["insertion"] = "sequence_not_extractable_in_test"
    mn, mx = kw.get("range", (0, 0))
    truth = {"insertion": ins, "upstream": up, "downstream": dn,
             "reference_empty": ref[refw_s:refw_e],
             "reconstructed_empty": up + dn}
    return {
        "_truth": truth,
        "schema_version": "1.0", "occurrence_id": occ, "event_id": "E_" + occ,
        "record_revision": 1,
        "carrier": {"source_id": carrier_src, "original_insert_start0": s,
                    "original_insert_end0": e, "insert_start0": s,
                    "insert_end0": e,
                    "canonical_insert_strand_in_carrier":
                        kw.get("canon_strand", "+"),
                    "interval_definition": "catalogue_excised_allele"},
        "reference": {"source_id": B, "junction0": j,
                      "reference_to_carrier_strand": kw.get("ref_strand", "+"),
                      "selection_evidence_id": "T"},
        "breakpoint": {"status": "unique" if (mn, mx) == (0, 0)
                       else "equivalent_range",
                       "min_joint_shift_bp": mn, "max_joint_shift_bp": mx,
                       "chosen_joint_shift_bp": 0,
                       "equivalence_evidence_id": "T"},
        "sequences": seqs,
        "quality": {"official": True, "analysis_set": "conservative",
                    "localization_status": "t",
                    "reconstructed_sequence_status": "t",
                    "candidate_search_truncated": False,
                    "core_recoverable": True, "delta_band": "t",
                    "event_delta_bp": 0, "auditA_verdict": "t",
                    "exclusion_reasons": []},
        "missing_sequence_reasons": miss,
    }


CASES = [
    ("ok_baseline", {}, "PASS"),
    ("reverse_canonical_insert", {"canon_strand": "-"}, "PASS"),
    ("reverse_reference_alignment", {"ref_strand": "-"}, "PASS"),
    ("real_snp_vs_reference", {}, "PASS"),
    ("window_offset_wrong", {"bad_junction_offset": True}, "FAIL"),
    ("wrong_assembly_same_contig",
     {"carrier_src": "SRC|ASM_B|contig1"}, "FAIL"),
    ("coordinates_off_by_one", {"ins_start": S0 + 1}, "FAIL"),
    ("join_order_reversed", {"reverse_join": True}, "FAIL"),
    ("invented_native_locus", {"native_locus": "ASM_A:contig1:0-100"}, "FAIL"),
    ("non_equivalent_range", {"range": (0, 5)}, "FAIL"),
    ("missing_sequence", {"drop_insertion": True}, "NOT_VERIFIABLE"),
]


def main():
    root = tempfile.mkdtemp(prefix="dbval_")
    try:
        car, ref = build(root)
        recs = [record("O%02d" % i, car, ref, **kw)
                for i, (_, kw, _) in enumerate(CASES)]
        with open(root + "/db/occurrences.jsonl", "w") as fh:
            for r in recs:
                fh.write(json.dumps({k: v for k, v in r.items()
                                     if k != "_truth"}) + "\n")
        # THE STORED BASES ARE ALWAYS THE TRUE ONES. A case mutates what the
        # record CLAIMS (coordinates, source, join order, offsets) and leaves
        # the sequence honest, so a FAIL means the validator noticed the claim
        # disagrees with the corpus -- not that a doctored FASTA disagrees
        # with itself. Generating the FASTA from the declared coordinates
        # instead would make every coordinate lie self-consistent and the
        # whole suite would pass vacuously.
        with open(root + "/db/sequences.fna", "w") as fh:
            for r in recs:
                for role, o in r["sequences"].items():
                    if o:
                        fh.write(">%s\n%s\n"
                                 % (o["sequence_id"], r["_truth"][role]))
        with open(root + "/v2/occurrences_flagged.tsv", "w") as fh:
            fh.write("occurrence_id\tofficial_fasta\tanalysis_set\n")
            for r in recs:
                fh.write("%s\tyes\tconservative\n" % r["occurrence_id"])

        env = dict(os.environ)
        for phase in ("records", "refrecords"):
            cmd = [sys.executable, VAL, "--pass", phase, "--shard", "0",
                   "--nshard", "1", "--out", root + "/work", "--w", root,
                   "--db", root + "/db", "--v2", root + "/v2",
                   "--src", root + "/src", "--replay-flank", "20"]
            p = subprocess.run(cmd, capture_output=True, text=True, env=env)
            if p.returncode != 0:
                print(p.stdout, p.stderr)
                raise SystemExit("validator crashed in pass " + phase)

        status = {}
        for fn in sorted(os.listdir(root + "/work")):
            if not (fn.startswith("val_") and fn.endswith(".tsv")):
                continue
            with open(root + "/work/" + fn) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    status.setdefault(f[0], []).append(
                        (f[1], f[2], f[3], f[4]))

        # also exercise the source-level duplicate-header check
        ps = subprocess.run(
            [sys.executable, VAL, "--pass", "sources", "--shard", "0",
             "--nshard", "1", "--out", root + "/work", "--w", root,
             "--db", root + "/db", "--v2", root + "/v2", "--src", root + "/src"],
            capture_output=True, text=True, env=env)
        dup_ok = False
        with open(root + "/work/src_000.tsv") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0].endswith("ASM_C|contig1") and f[2] == "FAIL" \
                        and "duplicate_fasta_record_id" in f[3]:
                    dup_ok = True

        npass = 0
        for i, (name, _, want) in enumerate(CASES):
            occ = "O%02d" % i
            rows = status.get(occ, [])
            has_fail = any(r[2] == "FAIL" for r in rows)
            has_nv = any(r[2] == "NOT_VERIFIABLE" for r in rows)
            if want == "PASS":
                ok = rows and not has_fail and not has_nv
            elif want == "FAIL":
                ok = has_fail
            else:
                ok = has_nv and not has_fail
            npass += ok
            why = [r for r in rows if r[2] != "PASS"]
            print("%-32s want=%-15s %s   %s"
                  % (name, want, "ok" if ok else "WRONG",
                     "" if ok else why[:3] or rows[:3]))
        print("%-32s want=%-15s %s"
              % ("duplicate_fasta_record_id", "FAIL",
                 "ok" if dup_ok else "WRONG"))
        npass += dup_ok
        total = len(CASES) + 1
        print("\n%d/%d negative-control cases behaved as specified"
              % (npass, total))
        if npass != total:
            sys.exit(1)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
