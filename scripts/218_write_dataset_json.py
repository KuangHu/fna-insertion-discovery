#!/usr/bin/env python3
"""dataset.json: the version, the counts, and the decisions this export made.

Counts are COMPUTED from the files just written, never copied from a previous
report. The spec's migration baseline is recorded separately, as a baseline to
compare against -- not as a fact asserted by this file.
"""
import hashlib
import json
import os
import subprocess
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
V2 = W + "/inferred_empty_v2"
OUT = W + "/insertion_db_v1"
FILES = ["source_manifest.jsonl", "occurrences.jsonl", "sequences.fna",
         "occurrences_view.tsv", "occurrences_browse.tsv",
         "reference_evidence.jsonl", "equivalence_evidence.jsonl"]
BASELINE = {"input_records": 177503, "input_events": 109709,
            "official_occurrences": 140843, "official_events": 84114,
            "conservative_occurrences": 130330, "conservative_events": 79105,
            "review_occurrences": 10513, "non_official_records": 36660}


def fsha(p, buf=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def git(*a):
    try:
        return subprocess.check_output(
            ["git", "-C", os.path.dirname(os.path.dirname(
                os.path.abspath(__file__)))] + list(a),
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


def main():
    counts = {"records": 0, "events": set(), "official": 0,
              "official_events": set(), "conservative": 0,
              "conservative_events": set(), "review": 0, "withheld": 0,
              "sequence_objects_present": 0, "sequence_objects_null": 0}
    with open(OUT + "/occurrences.jsonl") as fh:
        for line in fh:
            d = json.loads(line)
            counts["records"] += 1
            ev = d["event_id"]
            counts["events"].add(ev)
            q = d["quality"]
            if q["official"]:
                counts["official"] += 1
                counts["official_events"].add(ev)
            a = q["analysis_set"]
            if a == "conservative":
                counts["conservative"] += 1
                counts["conservative_events"].add(ev)
            elif a == "flagged_for_review":
                counts["review"] += 1
            else:
                counts["withheld"] += 1
            for r in ("insertion", "upstream", "downstream",
                      "reference_empty", "reconstructed_empty"):
                if d["sequences"].get(r):
                    counts["sequence_objects_present"] += 1
                else:
                    counts["sequence_objects_null"] += 1
    nseq = 0
    with open(OUT + "/sequences.fna") as fh:
        for line in fh:
            if line[0] == ">":
                nseq += 1
    nsrc = sum(1 for _ in open(OUT + "/source_manifest.jsonl"))

    actual = {k: (len(v) if isinstance(v, set) else v)
              for k, v in counts.items()}
    actual["sequence_records_in_fasta"] = nseq
    actual["source_contigs"] = nsrc

    ds = {
        "schema_version": "1.0",
        "dataset_id": "fna_insertion_occurrence_db_v1",
        "built_utc": subprocess.check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "synthetic": False,
        "coordinate_convention": "0-based-half-open",
        "stored_orientation": "native-contig-forward",
        "sequence_normalization":
            "uppercase_DNA_IUPAC_no_whitespace_no_wildcards",
        "files": {f: {"sha256": fsha(OUT + "/" + f),
                      "bytes": os.path.getsize(OUT + "/" + f)}
                  for f in FILES if os.path.exists(OUT + "/" + f)},
        "inputs": {
            "occurrences_flagged_tsv": {
                "path": os.path.relpath(V2 + "/occurrences_flagged.tsv", W),
                "sha256": fsha(V2 + "/occurrences_flagged.tsv")},
            "inferred_empty_fna": {
                "path": os.path.relpath(V2 + "/inferred_empty.fna", W),
                "sha256": fsha(V2 + "/inferred_empty.fna")},
        },
        "code": {"git_commit": git("rev-parse", "HEAD"),
                 "git_dirty": bool(git("status", "--porcelain")),
                 "scripts": ["214_export_manifest.py", "215_export_records.py",
                             "216_validate_db.py", "217_synth_db_validation.py"]},
        "parameters": {
            "carrier_flank_nominal_bp": 2000,
            "carrier_flank_clipped_at_contig_ends": True,
            "reference_window_rule":
                "mirror the carrier's own up/down flank lengths about the "
                "reference junction, clipped to the reference contig",
            "reconstruction_method": "concat_carrier_flanks_v2",
            "chosen_joint_shift_bp": 0,
        },
        "decisions": [
            "chosen_joint_shift_bp is 0 for every record: v2's chosen_shift is "
            "a REPORTED equivalent shift, not an applied one. Verified by "
            "replaying 3,012 official rows against the original FNA. The "
            "nonzero value is preserved as v2_reported_equivalent_shift_bp.",
            "sequence_accession_version is null: the pipeline holds assembly "
            "accessions, not per-contig accessions. Not fabricated.",
            "topology is unknown for every source: circularity was never "
            "determined. Circular-spanning occurrences are therefore out of "
            "scope for v1 rather than silently accepted.",
            "occurrences_view.tsv carries the inline sequence columns the spec "
            "lists; occurrences_browse.tsv is the same table minus those five "
            "columns, for routine viewing. The browse table is a strict column "
            "subset, not a second source of truth.",
        ],
        "counts_measured_from_this_build": actual,
        "migration_baseline_for_comparison": BASELINE,
        "limitations":
            "Source validation proves the stored bases are the bases at the "
            "declared coordinates of the declared assembly. It is NOT a "
            "measure of insertion accuracy and the two are never combined.",
    }
    with open(OUT + "/dataset.json", "w") as fh:
        json.dump(ds, fh, indent=2)
    print(json.dumps({"counts": actual, "baseline": BASELINE}, indent=2))
    diffs = [(k, actual.get(m), v) for k, m, v in (
        ("input_records", "records", BASELINE["input_records"]),
        ("input_events", "events", BASELINE["input_events"]),
        ("official_occurrences", "official", BASELINE["official_occurrences"]),
        ("official_events", "official_events", BASELINE["official_events"]),
        ("conservative_occurrences", "conservative",
         BASELINE["conservative_occurrences"]),
        ("conservative_events", "conservative_events",
         BASELINE["conservative_events"]),
        ("review_occurrences", "review", BASELINE["review_occurrences"]),
        ("non_official_records", "withheld",
         BASELINE["non_official_records"])) if actual.get(m) != v]
    if diffs:
        print("\nBASELINE DISAGREEMENT (measured vs spec baseline):")
        for k, a, b in diffs:
            print("   %-28s measured %s  baseline %s" % (k, a, b))
        sys.exit(1)
    print("\nall measured counts match the migration baseline")


if __name__ == "__main__":
    main()
