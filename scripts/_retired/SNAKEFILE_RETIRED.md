# The Snakefile is retired — it could not run

`Snakefile` (now `Snakefile.legacy`) described the ORIGINAL stage-10/40/50
pipeline. The SIMPLIFICATION work moved 13 scripts to `scripts/_retired/` and
the Snakefile was never updated with them, so on a fresh checkout it was broken
in two independent ways:

| | |
|---|---|
| `rule cluster` | calls `scripts/10_skani_cluster.py` — in `_retired/` |
| `rule merge`   | calls `scripts/40_merge_catalog.py` — in `_retired/` |
| `rule annotate`| requires `stage40/mmseqs_rep_seq.fasta`, which **no rule produces** — the DAG cannot be built at all |

It is kept for reference and is not an entry point. The pipeline that actually
runs is in `result_note/01_fna_to_insertion_call.md` §"The whole path";
every stage is submitted through `slurm/analysis.sh`.
