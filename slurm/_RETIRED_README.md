# `_retired_download_species.sh` — DO NOT RESUBMIT

Retired 2026-09-11. Downloads must not go through SLURM.

`cf1` is `OverSubscribe=EXCLUSIVE`, so this script billed a **whole 64-core
node** to perform pure network I/O: 225 CPU-hours for K. pneumoniae, 40 for
A. baumannii, with ~0 CPU actually used.

Use `tools/download_on_transfer_node.sh` on `lrc-xfer.lbl.gov` instead. The
five-species queue (9,863 genomes, 35 GB) ran there at **0 SU with 0 failures**.

See `docs/SIMPLIFICATION.md` §16.
