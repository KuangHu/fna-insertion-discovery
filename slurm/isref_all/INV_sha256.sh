#!/bin/bash
#SBATCH --job-name=isall_INVsha
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/global/scratch/users/kh36969/fna_ins_discovery/isref_all/inventory/%x_%j.out
# sha256 of the partial files left by the 27725308 tasks that SLURM cancelled
# at the time limit (shards 50, 52) and of shard 55's two hits files (gzip hit
# an I/O error, the done marker was written anyway).  Recorded in
# docs/legacy_shard_inventory.md.
#
# READ ONLY: the legacy directories are chmod -R a-w; this job writes only
# $W/inventory/.  Submit from the repo root:
#   mkdir -p $W/inventory
#   sbatch --export=ALL,GIT_COMMIT=$(git rev-parse HEAD) slurm/isref_all/INV_sha256.sh
set -euo pipefail
: "${GIT_COMMIT:?GIT_COMMIT not set}"
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
OUT="$W/inventory/sha256_${SLURM_JOB_ID:?run under SLURM}.txt"
FILES=(
  copies/shard_50.tsv
  hits/shard_50.raw
  hits/shard_50.tsv
  hits/shard_52.raw
  hits/shard_55.tsv.gz
  hits/shard_55.tsv
)
cd "$W"
{
  echo "# git_commit_submit=$GIT_COMMIT host=$(hostname) job=$SLURM_JOB_ID start=$(date -Is)"
  for f in "${FILES[@]}"; do echo "# before: $(stat -c '%s %Y %A' "$f") $f"; done
} > "$OUT.tmp"
printf '%s\0' "${FILES[@]}" | xargs -0 -P 6 -n 1 sha256sum > "$OUT.sums"
[ "$(wc -l < "$OUT.sums")" -eq "${#FILES[@]}" ] || { echo "FATAL: expected ${#FILES[@]} sums"; exit 1; }
{
  sort -k2 "$OUT.sums"
  for f in "${FILES[@]}"; do echo "# after: $(stat -c '%s %Y %A' "$f") $f"; done
  echo "# end=$(date -Is)"
} >> "$OUT.tmp"
rm -f "$OUT.sums"
mv "$OUT.tmp" "$OUT"
cat "$OUT"
