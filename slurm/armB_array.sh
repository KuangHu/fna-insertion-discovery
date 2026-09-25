#!/bin/bash
#SBATCH --job-name=armB_graph
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-49
#SBATCH --output=%x_%A_%a.out
# One array task per ANI cluster: graph build -> bubbles -> genotype -> refine.
#
#   sbatch --export=CLUSTERDIR=/path/stage10,OUT=/path/armB slurm/armB_array.sh
set -euo pipefail
: "${CLUSTERDIR:?set CLUSTERDIR}" ; : "${OUT:?set OUT}"
REPO=/global/home/users/kh36969/fna_based_mgefinder_project
CORE=/global/home/users/kh36969/.conda/envs/fnains
SV=/global/home/users/kh36969/.conda/envs/fnains_sv

mapfile -t MANIFESTS < <(ls "$CLUSTERDIR"/clusters/*.manifest | sort)
N=$SLURM_ARRAY_TASK_COUNT ; I=$SLURM_ARRAY_TASK_ID
for ((k=I; k<${#MANIFESTS[@]}; k+=N)); do
  m="${MANIFESTS[$k]}" ; cid=$(basename "$m" .manifest)
  od="$OUT/$cid" ; mkdir -p "$od"
  [ -s "$od/armB_events.tsv" ] && continue

  PATH="$CORE/bin:$PATH" python3 "$REPO/scripts/30_armB_graph.py" \
      --manifest "$m" --outdir "$od" --cluster-id "$cid" \
      --threads "${SLURM_CPUS_PER_TASK:-16}" --reuse \
      > "$od/stage30.log" 2>&1 || { echo "FAIL30 $cid" >> "$OUT/failures.txt"; continue; }

  PATH="$SV/bin:$PATH" python3 "$REPO/scripts/31_armB_refine.py" \
      --pairs "$od/armB_pairs.tsv" --clusters "$CLUSTERDIR/clusters.tsv" \
      --outdir "$od" --threads "${SLURM_CPUS_PER_TASK:-16}" \
      > "$od/stage31.log" 2>&1 || echo "FAIL31 $cid" >> "$OUT/failures.txt"
done
echo "task $I done"
