# Allele verification test module

The reusable module is `tools/classify_allele_evidence.py`. It consumes
anchor-pair mappings and original genomes, extracts both alleles, and writes
evidence and boundary status separately. It does not perform anchor discovery;
the local reference mapper is `tools/verify_filled_empty.py`.

Requires Python 3 and Biopython (tested with 1.88). Custom data:

```text
python tools/classify_allele_evidence.py --manifest genomes.tsv --verification verification.tsv --outdir allele_results
```

The tab-separated manifest has `genome_id` and `fasta_path` columns. Relative
FASTA paths resolve against the manifest directory. Genome and contig IDs must
match the verification table. Duplicate genome/contig IDs are rejected.

Required verification columns:

| Columns | Meaning |
|---|---|
| id, copy_id | Unique filename-safe row ID; biological copy label |
| source_genome, source_contig | Source sequence identifiers |
| source_start, source_end, source_insert_bp | Final copy interval and its length |
| boundary_refined | 0 or 1 from the copy table |
| target_genome, target_contig | Comparison sequence identifiers |
| target_start, target_end, target_strand | Inward anchor endpoints; strand relative to source |
| observed_gap_bp | target_end minus target_start; may be negative for overlapping anchors |
| qualifying_pairs | Number of qualifying mapped anchor pairs |
| status | Upstream mapping category; `no_comparator` means no comparison genome |

Coordinates are zero-based and half-open. The source is always extracted in
genomic orientation; reverse comparison windows are reverse-complemented before
alignment. Target start/end describe the genomic endpoints even on reverse
strand. Negative gaps do not define a biological negative-length allele.
Unmapped rows may use `.` for unused target fields. Multiple comparator rows
for the same copy need different row IDs to avoid overwriting evidence files.

Outputs include all input rows in `uniform_evidence.tsv`, observed nonempty
middle sequences in `uniform_middle_sequences.fna`, per-pair alignments,
`uniform_summary.json`, and `uniform_report.html`. Zero-length middles are
explicit in the table; undefined/overlapping middles are `.`. No original
candidate coordinates are changed. Full input validation precedes output writes.

The default no-argument invocation reruns the existing local dataset. The Python
API `run(rows, genomes, outdir)` accepts dictionaries of genome/contig sequences;
`compare(source_window, oriented_comparison_window, start, end)` returns metrics
and an alignment. Call coordinates for `compare` are window-relative.

## Automated tests

```text
python -m unittest discover -s tests -p "test_allele*.py" -v
python -m unittest discover -s tests -v
```

Tests use deterministic synthetic sequences and temporary directories, with no
downloads, minimap2, or local reference files required. They cover short sequence
retention, reverse genomic comparisons, zero and overlapping anchor gaps,
identical alleles, gaps outside the call, unresolved/multiple placements, missing
comparators, invalid inputs, and an end-to-end custom-manifest CLI run. The
separate existing POSIX mock-aligner test is skipped on Windows.

## Interpretation limits

Support uses the documented uniform screening rule: one source-only gap covers
at least 90% of the called interval, aligned-base identity is at least 95%, and
both source flanks contribute at least 200 aligned bases. These are operational
thresholds, not calibrated probabilities. Repeated junctions can have equivalent
optimal alignments. A short sequence change can split a large gap, failing the
single-gap threshold despite a real long/short difference; a regression test
preserves this limitation. `no_dominant_long_short_support` never means proven
absence. No motif, TSD, replacement mechanism, or insertion polarity is assigned.
