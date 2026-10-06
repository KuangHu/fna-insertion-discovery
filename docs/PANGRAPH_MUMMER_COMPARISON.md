# Local Pangraph and independent MUMmer comparison

Inputs: the same original R6 and TIGR4 complete references used for Arm B.
Production Arm B and locked Layer 2 were not modified. The retired Pangraph
adapter was not revived: the current 1.4.0 node/path JSON is parsed explicitly.

## Reproduce

In WSL, run `python3 tools/local_pangraph_mummer.py`. It requires the official
Pangraph 1.4.0 binary in `env/bin/pangraph`, and MUMmer 3.23 (`nucmer`,
`delta-filter`, `show-coords`, `show-diff`) on PATH. Then run
`python tools/compare_pangraph_mummer.py`. Original graphs, command arguments,
logs, raw delta and coordinates, and one-to-one filtered mappings are retained
under `local_debug/pangraph_mummer/`.

Pangraph used `build --verify --circular -j 2` with default alignment parameters
and both input orders. Both builds passed exact genome reconstruction checks.
The default Pangraph backend uses minimap2 internally, so Pangraph is a different
graph construction method, not an entirely independent alignment engine.
MUMmer `nucmer --maxmatch` supplies the independent alignment evidence here;
SyRI was not run.

## Discovery comparison

Both graphs contain 368 blocks and 166 blocks occurring once in each genome.
Taking adjacent single-copy shared blocks in matching orientation gives 160
comparable inter-block intervals: 106 below 500 bp net difference, 43 from
500 to 5000 bp, and 11 over 5000 bp. Coordinates are identical across input
orders. Five order/orientation-incompatible anchor pairs and the origin-crossing
pair were excluded, not called negative. Within-block variants were not called
separately (maximum observed within-block deletion 93 bp).

A context link requires intervals to overlap or lie within 200 bp in BOTH
genomes. By this criterion 44/46 Arm B primary candidates have a Pangraph
counterpart, sometimes in a larger interval. All 43 Pangraph primary intervals
have an Arm B context link. This is many-to-many regional concordance, not a
recall estimate or proof that events and their boundaries are identical.
Arm B B000003 and B000013 lack a link under this extraction rule; this does not
establish that the full Pangraph graph lacks the relevant sequence.

## Four disputed loci: independent alignment

MUMmer evidence comes from all alignments, not only a one-to-one filter. For
each known locus in both genomes, local forward flanking blocks with >=200 bp
alignment and >=95% identity were sought within 2 kb. Candidate pairs were
ranked by distance to the known endpoints. These are alignment extents, not
biochemical boundaries; input candidates constrain the search location.

| Arm B candidate | Pangraph net difference | MUMmer evidence |
|---|---:|---|
| B000002 | 912 bp | overlapping local repeat alignments; no nonoverlapping flanking pair under this rule, remains unresolved |
| B000012 | 1172 bp | flanking blocks bracket 999 bp in R6 versus 2171 bp in TIGR4: net 1172 bp |
| B000019 | 832 bp | broader interval is 1174 bp in R6 versus 1989 bp in TIGR4: net 815 bp, not an exact boundary confirmation |
| B000036 | 1719 bp | flanking blocks bracket 1 bp in R6 versus 1709 bp in TIGR4: net 1708 bp |

For B000002, R6 intervals 75238–82431 and 81868–92717 align to overlapping
TIGR4 intervals 81348–88543 and 87068–97909 (0-based half-open). This is
consistent with the placement ambiguity already flagged by Layer 2; it is not
a reason to force a unique insertion boundary. One-to-one filtered MUMmer
output also contains remote repeat matches with apparent inversion labels,
which are not accepted here as proven rearrangements.

Conclusion: Pangraph provides a useful stable supplementary representation on
this two-genome test, but does not establish superior discovery sensitivity.
Independent MUMmer evidence strengthens the 1708 bp difference at B000036 and
supports a complex length difference at B000012; B000019 retains boundary
uncertainty and B000002 retains repeat-placement ambiguity. Keep Arm B →
Layer 2 as the main route and retain these cross-checks as separate evidence.

Outputs: `report.html`, `summary.json`, both `*_intervals.tsv`,
`candidate_links.tsv`, `disputed_mummer.tsv`, and per-locus regional alignment
tables. Raw sequences/graphs/alignments remain alongside these reports.
