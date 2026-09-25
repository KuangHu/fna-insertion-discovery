---
name: long-jobs-need-slurm
description: On this cluster, anything longer than a few minutes must go through sbatch — background shells die with the session
metadata:
  type: feedback
---

Run long jobs with `sbatch`, never as a background shell from the session. Backgrounded commands here are tied to the session's process lifetime: `setsid nohup ... &` did NOT survive, and neither did already-running children.

Burned on 2026-08-27: an ISEScan run over 20 genomes was launched as a background shell on a compute node, ran ~75 minutes, and was killed at session teardown with **zero** `.tsv` produced. All 20 genomes had started; all 20 were lost.

**Why:** the harness reaps the process group, so a multi-hour job silently loses everything with no partial output.
**How to apply:** if a job might exceed a few minutes, write a `slurm/*_array.sh` in the repo convention (see `slurm/armA_array.sh`: `--account=pc_rubinlab --partition=lr6 --qos=lr_normal`, `FNA_LIST`/`OUT` via `--export`, strided `awk 'NR%n==i'`, and a skip-if-output-exists guard so requeues resume). Then poll `squeue`. Live example: `slurm/isescan_array.sh`, submitted as job 25306752 over 20 genomes with `--export=FNA_LIST=...,OUT=/global/scratch/users/kh36969/fna_ins_discovery/isescan_run`, normalising results to `isescan_run/tsv/<sample>.tsv`. It also gives each task its own CWD, because ISEScan drops scratch `<dir>_<name>.list` files into the working directory and concurrent tasks otherwise collide. Note ISEScan specifically drives `phmmer --max`, which is very slow (tens of min/genome) and barely scales with `--cpu` — parallelise across genomes, not threads.

See [[fna-ins-discovery-state]].
