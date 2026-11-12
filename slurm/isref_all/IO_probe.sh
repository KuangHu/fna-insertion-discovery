#!/bin/bash
#SBATCH --job-name=isall_IOprobe
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=00:50:00
#SBATCH --exclude=n0052.cf1
#SBATCH --output=/global/scratch/users/kh36969/fna_ins_discovery/isref_all/ioprobe/logs/%x_%j.out
# Read-only I/O probe of the genome FASTA files (diagnostic; user-approved
# 2026-11-11 after V0 job 27788052 hit 25 I/O errors on n0059.cf1).
# Writes only $W/ioprobe/<job>/.  Submit from the repo root:
#   mkdir -p $W/ioprobe/logs
#   sbatch --export=ALL,GIT_COMMIT=$(git rev-parse HEAD) slurm/isref_all/IO_probe.sh
set -euo pipefail
: "${GIT_COMMIT:?GIT_COMMIT not set}"
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
REPO=/global/home/users/kh36969/fna_based_mgefinder_project
OUT="$W/ioprobe/${SLURM_JOB_ID:?run under SLURM}"
echo "# git_commit_submit=$GIT_COMMIT host=$(hostname) job=$SLURM_JOB_ID start=$(date -Is)"
echo "# kernel $(uname -r); lustre client: $(lfs --version 2>&1 | head -1)"
findmnt -T /global/scratch -o TARGET,FSTYPE,OPTIONS || true
echo "# lfs df at start"; timeout 120 lfs df -h /global/scratch || echo "lfs df failed rc=$?"
echo "# lfs check servers"; timeout 120 lfs check servers 2>&1 | grep -v ' active\.$' | head -50 || true
"$PY" -I "$REPO/scripts/249_io_probe.py" --corpus "$V/corpus_v2.tsv" --out "$OUT"
echo "# lfs check servers at end"; timeout 120 lfs check servers 2>&1 | grep -v ' active\.$' | head -50 || true
echo "# end=$(date -Is)"
