---
name: fna-ins-discovery-paths
description: Where the FNA-only mobile-element pipeline's code, run outputs, envs and test genomes live
metadata:
  type: reference
---

Code: `/global/home/users/kh36969/fna_based_mgefinder_project` (NOT a git repo — no history to recover from).
Run outputs / workdir: `/global/scratch/users/kh36969/fna_ins_discovery` (`workdir:` in `config/config.yaml`).
Test genomes: `/global/scratch/users/kh36969/cross_ref_is_runs/results/escherichia_coli/_ncbi_symlinks/*.fna` (200 staged E. coli).

Runtime PATH that makes stages 10/20/30/30b/40 work with no conda solve:
`export PATH=$PWD/env/bin:$HOME/.conda/envs/claude-env/bin:$PATH`
- `env/bin/` holds static binaries: minigraph 0.21-r606, gfatools 0.5-r296, skani 0.3.2, pangraph 1.4.0.
- `claude-env` supplies minimap2, mmseqs, samtools, nucmer, prodigal, hmmsearch.
- The `fnains` core env from `setup_env.sh` was never built and is not needed.
- Built and working: `fnains_sv` (svim-asm), `fnains_annot` (isescan.py, barrnap, hmmsearch, seqkit), `fnains_util` (nucdiff, seqkit).

See [[fna-ins-discovery-state]].
