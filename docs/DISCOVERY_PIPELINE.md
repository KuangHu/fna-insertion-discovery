# Arm B → Layer 2 workflow

The supported primary route is structural candidate discovery by Arm B,
followed by nucleotide allele reconstruction with the existing locked Layer 2.
Arm A provides complementary multicopy evidence. Pangraph and MUMmer provide
optional supplementary evidence and never remove primary candidates.

## Run

Requires Python 3 with Biopython, minigraph, gfatools, minimap2 and samtools.
Use an existing homologous assembly panel, with one uncompressed FASTA filename
per manifest line. Relative paths resolve against the manifest directory.
The runner does not perform ANI grouping. Choose a new output directory for
each run to prevent stale results from being reused.

`GIT_COMMIT` must be set to the commit being run (G-4); both the runner and
`30_armB_graph.py` refuse to start without it.

```bash
GIT_COMMIT=$(git rev-parse HEAD) python3 tools/discovery_pipeline.py \
  --manifest panel.manifest --outdir results/panel --panel-id panel \
  --minigraph env/bin/minigraph --gfatools env/bin/gfatools \
  --minimap2 env/bin/minimap2 --threads 2
```

For two complete single-contig assemblies with distinct FASTA record IDs, add
`--supplement --pangraph env/bin/pangraph`. This runs Pangraph with exact input
reconstruction verification and independently aligns the two genomes with
MUMmer (`nucmer --maxmatch`, raw and one-to-one-filtered results retained).
Add `--circular` only for circular genomes. Supplementary comparison currently
supports two genomes, whereas the primary route supports multi-genome panels.
MUMmer executables are required when requesting supplementary evidence.

## Outputs and interpretation

- `armB/`: original graph, bubbles, candidate tables, sequences, and parked
  events larger than 5 kb. Arm B default scope is a 500–5000 bp length delta.
  Arm B's contract is SPEC item M1-B (2026-11-11):
  - its outputs are usable only if `armB/armB_complete.json` exists with
    `status: complete`. The marker is written last and records counts, output
    sha256 and git provenance; a failed run leaves no final output files;
  - with `--reuse`, a cached graph or call file is used only if it was made
    from the same inputs (and, for calls, the same graph) by the same
    minigraph (measured version and executable sha256) with the same
    arguments. The marker's `cache` field says, for the graph and every call
    file, whether it was reused and carries the generator and creation record
    (time, commit, host) of the run that made it;
  - each genome's `.fai` is checked line by line against its FASTA before
    any base is read (`scripts/lib/fasta.py`); an index that does not describe
    the file is fatal;
  - every genome at a bubble is counted as empty, filled, intermediate or
    missing (`n_*` columns); bubbles below the size gate are counted in the
    marker, and bubbles with more than two length classes are flagged
    `multiallelic`;
  - the insert is the difference between the empty and filled reference
    genomes' own graph walks (`empty_ref_genome`, `filled_ref_genome`). If its
    length differs from `insert_size` by more than `--len-tol`, the event is
    flagged `insert_len_mismatch`;
  - `junction_L` / `junction_R` are 0-based half-open contig coordinates of
    `empty_ref_genome`. They are given only when `frame_status` is `verified`,
    meaning that genome's own allele, with 50 bp of source/sink context, was
    found exactly once in its whole contig (both strands), on the call's
    strand, overlapping its minigraph call interval, and the whole call
    interval matches the graph path at that position. Otherwise they are -1
    and `frame_status` gives the reason: `unverified` (no exact copy),
    `ambiguous` (more than one), `off_call_interval` (the only exact copy is
    elsewhere) or `interval_mismatch` (the interval's own sequence does not
    line up). The raw call interval is always given in
    `empty_ref_call_start/end`. It is approximate: slack of up to ~570 bp was
    measured;
  - `parent_large_event_id` has been removed (it was never populated).
- Layer 2 seeds: only a `verified` frame is seeded, at junction ± 200 bp
  (`seed_status: seeded`). Every other frame has no verified junction and is
  kept as an unseedable candidate (`frame_<status>`), as are
  `genome_or_contig_not_in_manifest` and `seed_padding_outside_contig`. On
  one real E. coli panel (GCA_013745955.1_k24, 1,503 calls), 70.0% of calls
  verified, 27.3% had no exact copy (mostly a SNP near the site), 1.6%
  `interval_mismatch`, 1.0% `off_call_interval` and 0.1% `ambiguous`. Per
  event the verified share is higher, because any empty genome that verifies
  can serve as the frame.
- `layer2_seeds.tsv`: seeds as described above.
- `layer2/`: original Layer 2 allele sequences and decomposition/placement data.
- `layer2_tmp/`: retained anchors and mapping evidence from the locked script.
- `candidates.tsv`: all primary Arm B candidates, with separate Layer 2 status
  and prefixed reconstruction fields. Unseedable, unresolved and ambiguous
  candidates remain visible. An output-free Layer 2 locus is unresolved.
- `supplementary/`: verified Pangraph JSON, conservative inter-block interval
  table, independent MUMmer delta, coordinates and difference classifications.
- `report.html`, `summary.json`, `provenance.json`, `logs/`: reviewable report,
  statuses (including seed status and Arm B counts), input/executable/script
  SHA256, exact commands, and git provenance (`git_commit_runtime`,
  `git_dirty`, and `git_commit_submit` from `$GIT_COMMIT`, recorded as absent
  when not set).

Graph delta (`insert_size`), exported long fragment (`insert_len_direct`) and
reconstructed fragment (`layer2_inserted_len`) are distinct quantities.
Short intervals are retained verbatim; padded allele lengths are not target
motif lengths. Layer 2 event classes describe sequence representations and do
not determine evolutionary polarity or certify biochemical reaction sites.

Pangraph comparison uses adjacent blocks that occur once in both genomes and
have matching local order/orientation. Origin-crossing and reordered pairs are
excluded. Within-block variation remains in the graph but is not independently
called by this interval extractor. A Pangraph interval is supplementary regional
evidence, not a proven insertion. Pangraph's default alignment backend uses
minimap2; MUMmer provides the independent aligner comparison.

Primary failures stop the workflow. Supplementary failures return failure while
preserving completed primary output and explicitly marking supplementary status.
The runner leaves the locked Layer 2 algorithm unchanged; the historical
251-event benchmark is not part of the local smoke test.

For the two-order local comparison, the case-study tools remain available:
`tools/audit_armB.py`, `tools/audit_armB_layer2.py`, and
`tools/compare_pangraph_mummer.py`. Their R6/TIGR4-specific reports are distinct
from the general workflow outputs.

## Local validation

`python -m unittest discover -s tests -v` runs regression tests. The POSIX mock
aligner test requires WSL. Use the real R6/TIGR4 local manifests to exercise all
tools end-to-end; results and binary dependencies are intentionally not tracked
in Git. Download source accessions and hashes are recorded by local_debug.py.

The refactored workflow was rerun with both R6/TIGR4 input orders, including all
supplementary steps. Each order recovers 46 primary candidates and 12 parked
large events. Layer 2 statuses are 40 resolved / 1 ambiguous / 5 unresolved
(R6-first), and 42 / 1 / 3 (TIGR4-first). Both Pangraph runs recover 43 primary
inter-block intervals. `tools/verify_refactored_local.py` verified that original
Arm B tables/sequences, Layer 2 tables/sequences, and Pangraph JSON are identical
to the pre-refactoring audited outputs; forward MUMmer coordinates also match.
The complete WSL regression suite passed 45 tests. The Layer 2 implementation
was unchanged, and its historical 251-event benchmark was not rerun.
