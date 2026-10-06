# Local debugging on this computer

Setup date: 2026-10-05. Checkout baseline: `9316500b42c6c8094b29c2dfbd918cfe2299e273`.

The repository is at `F:\science\chatgpt\RNA_guided_editor\fna-insertion-discovery`.
It runs through the existing Ubuntu WSL2 distribution. SSH authentication to
GitHub was unavailable, so the checkout uses HTTPS. Changes have not been pushed.

## Installed and exercised

- Ubuntu Python 3.14 (standard library only for this test).
- minimap2 2.30-r1287, official release binary in `env/bin/minimap2`.
- samtools 1.22.1 installed in Ubuntu, for FASTA indexing.
- Arm A discovery, script 181 alignment expansion, and the new local record audit.

The full graph/annotation pipeline is NOT installed or tested. In particular,
this setup does not yet cover Arm B, skani, minigraph, MMseqs2 or HMM annotation.

## Rerun from PowerShell

```powershell
wsl -d Ubuntu -- bash -lc 'cd /mnt/f/science/chatgpt/RNA_guided_editor/fna-insertion-discovery && python3 tools/local_debug.py'
wsl -d Ubuntu -- bash -lc 'cd /mnt/f/science/chatgpt/RNA_guided_editor/fna-insertion-discovery && python3 -m unittest discover -s tests -p "test_local*.py" -v'
```

The runner reuses downloaded genomes, overwrites only its derived local-debug
outputs, and records source URLs, SHA-256 hashes, commands and wall-clock times.
The ignored `local_debug/` directory contains data, logs, seeds and results.

## Downloaded examples and observed output

The count table below records the initial smoke test, before the Arm A boundary
fix. Current rerun counts and boundary verification are recorded at the end.

| Reference | Assembly | Bases | Arm A candidate representatives |
|---|---|---:|---:|
| S. pneumoniae R6 | GCF_000007045.1 | 2,038,615 | 5 |
| S. pneumoniae TIGR4 | GCF_000006885.1 | 2,160,842 | 7 |
| E. coli K-12 MG1655 | GCF_000005845.2 | 4,641,652 | 7 |

Sources are NCBI RefSeq genomic FNA downloads; exact URLs and file hashes are in
`local_debug/download_manifest.json`.

Nineteen Arm A representatives supply the smoke-test queries. Each receives an
artificial DEBUG cluster ID: these are NOT production CDS families and are not
deduplicated across genomes. This intentionally exposes duplicate-seed hits.

| Target collection | Passing hits | Candidate position groups | Multi-label groups | Exact flank classes |
|---|---:|---:|---:|---:|
| Two S. pneumoniae genomes | 119 | 70 | 49 | 66 |
| E. coli K-12 | 36 | 36 | 0 | 36 |

All 155 hits were same-species in this small dataset. This tests the expansion
code path but does NOT demonstrate cross-species transfer or host-range breadth.
No independent insertion counts, RNA guidance, biological family assignments,
or empty-allele evidence are established by this run.

## Fixes and tests

Script 181 now preserves all passing alignments instead of suppressing records
with the same `(cluster, genome, contig, start//50)` key. It records query bounds,
strand and the coordinate convention, creates the output directory and fails
loudly on an aligner error. A CLI regression supplies two hits in the same old
50-bp bin and verifies both survive with their orientations. The same test
checks that a failing aligner cannot silently look like a completed empty run.

The original script 186 is preserved for comparison. `scripts/local_record_objects.py`
is a separate local audit with explicit input/output paths and a different schema:

- Stable sorting and a full deterministic tie-break for the display representative.
- Greedy complete-link interval grouping across all active groups, so an intervening
  nested hit does not split two surrounding equivalent intervals.
- No transitive overlap chains that merge intervals failing the group criterion.
- Separate hit, seed and family-label counts; actual event IDs in `all_seeds`.
- Full SHA-256 flank keys rather than 12-base display prefixes as identifiers.
- Incomplete or ambiguous-base flanks excluded from exact-key classes.
- Every original hit retained in the companion `.hits.tsv`; the runner verifies
  equality of the complete input/output record multisets.
- Host-locus and pairing-variant inference explicitly marked `not_inferred`.

Nine regression tests pass. They include score ties/input shuffling, nested
intervals, overlap chains, prefix collisions, incomplete flanks, multiple seeds
in one family, opposite query orientations and invalid intervals.

## Reading the output

`local_debug/expansion/spneumoniae_objects.tsv` is the representative site-group
view. `spneumoniae_objects.hits.tsv` retains every contributing alignment and
its own flanks. Equivalent files are present for E. coli.

The representative is selected by identity, then query coverage, then a stable
record key. It is only a display choice. Its boundary and flank are not validated
element boundaries. Start/end ranges expose alternatives. Opposing query strands
can occupy one physical group and remain in the evidence table. Flank classes
remain in target-contig orientation and are not reverse-complement canonicalized.

The clustering is a reproducible heuristic, not proof of physical insertion
identity; nested or closely spaced biological events still need inspection.

The original eight-seed S. aureus-to-S. pneumoniae case cannot be reproduced from
the Git checkout alone: its HPC catalogue, hit table and target assembly are not
included. The public examples validate local execution and specific code fixes,
not that original biological case.

## Arm A boundary repair (2026-10-05)

`project_target_boundary` now walks the forward PAF `cg` CIGAR, accounting for
insertions and deletions. It rejects missing/malformed CIGAR, out-of-alignment
positions, deletion interiors and insertion junctions with multiple possible
query coordinates. The existing integer conversion of median coordinates is
performed once in the representative frame. If either end is unresolved, the
whole original copy interval is retained and `boundary_refined=0`; the family
also records the unresolved projection count. Unaligned copies are no longer
mistakenly marked as refined.

Both representative FASTA and copy exports now call `resolved_copy` for the
same effective coordinates. `rep_copy`, the representative sequence and its
MD5 therefore agree with the corresponding copy row. Raw family length median
and consensus length remain distinct statistics, not representative lengths.

Validation:

- 20 regression tests, including direct gap projections, ambiguous endpoints,
  retained nominal intervals and the previous local expansion/grouping tests.
- Three-genome rerun: 19 representatives and 117 copies checked against original
  genome slices, oriented flanks, site FASTA and representative MD5. All pass,
  including 44 reverse-strand copies.
- 29 copies retain nominal boundaries (R6: 8, TIGR4: 15, E. coli: 6).
- Independent R6 CIGAR audit: of the former 11 internally checkable boundaries,
  8 remain on accepted refined copies and all agree with CIGAR projection.
  Three are excluded because the other end of their copy is unresolved; they
  are not counted as corrected, validated projections.
- Current expansion results: S. pneumoniae 127 hits / 86 candidate groups;
  E. coli 36 / 36. These are consequences of changed query sequences and do
  not establish increased biological accuracy or independent insertion counts.

The previous Arm A outputs are retained in `local_debug/armA_before_boundary_fix`.
Machine-readable verification: `local_debug/export_verification.json` and
`local_debug/boundary_projection_audit.json`. Run the sequence export check with
`python3 tools/verify_armA_exports.py` from the repository inside WSL.

### Filled / short allele audit

Run `python3 tools/verify_filled_empty.py` inside WSL, then run
`python tools/audit_allele_pairs.py` with Biopython available. These are independent
audit tools; neither rewrites Arm A calls or the locked Layer 2 implementation.

The 117 source copies were re-extracted with 500 bp on each side. This confirms
source-coordinate consistency, not biological insertion boundaries. Of 77
S. pneumoniae copies compared between R6 and TIGR4, 3 had positive short anchor
gaps (5, 9, 10 bp), 7 had overlapping anchors, 22 had similarly sized middle
intervals, 1 had another interval length, and 44 had unresolved anchors under
the stated thresholds. No exact zero-gap anchor pair was observed. E. coli's
40 copies have no second E. coli genome in this local dataset.

Independent global sequence alignments of the 3 short and 7 overlapping pairs
all contain a dominant source-only gap covering more than 99.4% of the called
interval. Aligned-base identities are 96.79–99.90%; they exclude gap columns.
Gap endpoints differ from nominal call endpoints by up to 33 bp. These results
support long/short allele differences but do not uniquely resolve breakpoints.
The three positive short-gap cases all have `boundary_refined=0`. Their
call lengths / observed short intervals / alignment gap lengths are respectively
1448/5/1443, 1717/9/1708, and 1717/10/1707 bp. The seven negative anchor gaps
must not be interpreted as negative-length alleles or proven TSDs.

The alignment audit uses one optimal global alignment (match +2, mismatch −3,
gap open −12, extension −0.1). Equivalent gap placements in repeats remain
possible. This is assembly-level evidence: no evolutionary polarity, independent
event count, transposition mechanism, or RNA guidance is established. Missing
anchor pairs are unresolved, not absence; similarly sized intervals do not
establish sequence identity. Unique placement means one qualifying pair among
the reported mappings and thresholds, not proven genome-wide uniqueness.

Results in `local_debug/filled_empty/`: `verification.tsv` (all copies),
`pair_alignment_audit.tsv` (ten sequence comparisons), `sequence_audit.html`
(readable architecture and evidence links), and per-pair FASTA/full alignments.

### Uniform candidate evidence policy

Run `python tools/classify_allele_evidence.py` after the anchor audit (Biopython
required). Every copy uses the same final copy-table interval. Every uniquely
paired locus, regardless of anchor-gap category, undergoes the same whole-window
alignment, including reverse-oriented comparison loci. Unmapped and no-comparator
copies remain explicit rows. Comparison middle sequences are retained verbatim
in source orientation; overlapping anchors have an undefined middle interval.

`evidence_status` is separate from `boundary_status`. Operational support requires
one source-only gap covering >=90% of the call, >=95% identity among aligned base
columns, and >=200 aligned bases in each source flank. These thresholds are a
transparent screening rule, not a calibrated probability or proof of mechanism.
No motif/TSD/replacement interpretation is assigned. Boundaries are not rewritten.

Local rerun: all 117 rows retained; 33 unique pairs aligned with identical rules;
10 supported long/short differences, 23 without dominant long/short support,
44 unresolved, 40 with no comparator. The former ten examples are recovered
without selecting them by their short/overlapping anchor status. Outputs:
`uniform_report.html`, `uniform_evidence.tsv`, `uniform_middle_sequences.fna`,
`uniform_summary.json`, and per-pair `*_uniform_alignment.txt`.

### Arm B real-genome run and audit

Arm B now runs locally under WSL: `python3 tools/local_armB.py`. This builds
R6/TIGR4 panels in both backbone orders, with default stage-30 settings (minimum
length difference 500 bp, maximum 5000 bp, 500 bp flank QC, two threads).
Official dependencies were built under `env/src`; minigraph commit
`2f569ebe3071fcb242f1bf0eabbb917941b47239` (0.21-r606), gfatools commit
`39c54d5dfe4678d4cfb75b37bd27d4fb21f4de72`. This is a two-reference smoke test,
not a population benchmark or complete test of panel selection.

Run `python tools/audit_armB.py` with Biopython to audit every output against
original genomes. Both orders yield 145 bubbles, 46 primary structural-allele
candidates and 12 parked large events. All 58 exported sequences per order
pass the >=95% local sequence identity/coverage check against the long carrier;
coordinates and export MD5 also pass. Large events have sequence/coordinate
checks, but not full long-versus-short global alignment.

All 46 primary pairs per order have complete genome-window FASTAs and global
alignments. The uniform dominant-gap rule supports 27 (R6-first) and 28
(TIGR4-first); the others are not called false positives. Shared sequence within
a bubble, complex differences, and imperfect homology can fail that rule.

The 46 primary exported fragments match MD5-for-MD5 in the two runs, but only
42 pairs match both R6 genomic endpoints within 100 bp. The four remaining
same-sequence candidates have shifted call boundaries, not demonstrated losses.
Graph path length differences disagree with actual called genomic interval
length differences at 15/46 and 16/46 events respectively (maximum 475/461 bp).
This distinction is real in the minigraph call BED and is retained by the
parser: a zero-length graph walk may still span several genomic bases.

Example `r6_first.B000000`: graph delta 1443 bp, matching the previously audited
R6/TIGR4 locus. Example `r6_first.B000001`: graph delta 502 bp but exported long
fragment 1316 bp, with an 814 bp short graph path. `diff_alleles` trims only
exact common prefixes/suffixes; its output is not necessarily a pure added
element. Keep graph delta, genomic delta, and exported fragment length separate.

No production Arm B calls or algorithm were changed. Findings support using
Arm B for candidate discovery followed by independent nucleotide/anchor
reconstruction; its graph coordinates and TSD/architecture annotations are not
mechanistically verified. Audit files: `local_debug/armB/audit.html`, `audit.tsv`,
`audit_summary.json`, `order_pairs.tsv`, and both panels' `audit/` directories.

### Arm B → locked Layer 2 backbone-order verification

Run `python3 tools/local_armB_layer2.py` in WSL (Biopython 1.86 installed), then
`python tools/audit_armB_layer2.py` for the report. Each run uses all 46 primary
events, seeding in the designated short carrier at `[junction_L-200,
junction_R+200)`, as documented in PIPELINE_RECTIFICATION.md. Layer 2 uses its
default 500 bp anchors and 7 kb interval cap. Its source was not modified;
SHA256 is recorded in `local_debug/armB/layer2_run.json`. No claim is made that
the separate historical 251-event regression benchmark was rerun.

R6-first resolves 41/46 loci (40 unique, 1 ambiguous); TIGR4-first resolves
43/46 (42 unique, 1 ambiguous). All 41 common loci agree in sequence-class
labels, 38 in reconstructed fragment length and MD5, and 37 in both absolute
junction coordinates. Reverse-complement/repeat-slide equivalence does not
rescue the three discordant fragment pairs. All three are length polymorphisms,
not simple insertion representations. A fourth length-polymorphism locus has
the same fragment but a 19 bp difference in its right junction.

All 29 shared insertion-target-retained/replacement representations agree in
fragment and absolute junction coordinates; 28 have unique placement in both
runs. The short allele includes ordinary padding and is not a target motif.

Previously order-sensitive loci:

| Candidate | Layer 2 result |
|---|---|
| B000002 | 912 bp fragment and junctions agree, but placement remains ambiguous in both runs |
| B000012 | length polymorphism in both; fragment 2477 vs 2650 bp, right junction differs 173 bp |
| B000019 | unresolved R6-first; length polymorphism TIGR4-first |
| B000036 | unresolved R6-first; 1708 bp replacement representation TIGR4-first |

B000008, B000015 and B000039 remain unresolved in both. Missing loci are
explicitly retained in the 46-row comparison report, even though Layer 2's
native loci.tsv only emits loci placed in at least two genomes.

168 emitted genome/allele occurrences were checked for length/hash and exact
occurrence in their listed original carrier (either orientation). All inferred
fragments were independently re-sliced from the long allele using
`insert_start_in_long_bp` and verified against `inserted_md5`. These checks show
sequence consistency, not uniquely correct placement. The results demonstrate
partial order stability; they do not prove biochemical junction accuracy or
complete convergence. No production boundaries or candidate lists were changed.

Outputs: `local_debug/armB/layer2_comparison.html`, `layer2_comparison.tsv`,
`layer2_comparison.json`, and each panel's `layer2/` original output files.
