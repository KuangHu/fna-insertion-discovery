#!/bin/bash
#SBATCH --job-name=isescan
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-19
#SBATCH --output=%x_%A_%a.out
# ISEScan for the stage-60 contingency table (Arm A x ISEScan).
#
# This MUST be a batch job, not a background shell. ISEScan drives
# `phmmer --max`, which disables the heuristic filters and is orders of
# magnitude slower than default phmmer -- tens of minutes per genome, and it
# barely scales with --cpu (measured: ~0.3-0.75 of a core per process on an
# otherwise idle 40-core node). A session-scoped background run gets killed
# long before it finishes.
#
#   sbatch --export=FNA_LIST=/path/list.txt,OUT=/path/isescan_run \
#          slurm/isescan_array.sh
#
# FNA_LIST is one .fna path per line. Each array task takes a strided slice, so
# the array size and the genome count do not have to match.
set -euo pipefail
: "${FNA_LIST:?set FNA_LIST}" ; : "${OUT:?set OUT}"
export PATH="/global/home/users/kh36969/.conda/envs/fnains_annot/bin:$PATH"
mkdir -p "$OUT/logs"

N=$SLURM_ARRAY_TASK_COUNT ; I=$SLURM_ARRAY_TASK_ID
awk -v n="$N" -v i="$I" 'NR%n==i' "$FNA_LIST" | while read -r fna; do
  [ -s "$fna" ] || continue
  s=$(basename "$fna"); s=${s%.gz}; s=${s%.fna}; s=${s%.fa}; s=${s%.fasta}
  # ISEScan drops scratch files (`<dir>_<name>.list`) into the CWD and keys its
  # output tree on the input path, so each task gets its own sandbox. Sharing
  # one CWD across concurrent tasks is what littered the last run.
  work="$OUT/work/$s"
  [ -s "$OUT/tsv/$s.tsv" ] && continue
  mkdir -p "$work" "$OUT/tsv"
  cp -f "$fna" "$work/$s.fna"
  ( cd "$work" && isescan.py --seqfile "$work/$s.fna" --output "$work/out" \
        --nthread "${SLURM_CPUS_PER_TASK:-4}" ) \
      > "$OUT/logs/$s.log" 2>&1 || { echo "FAIL $s" >> "$OUT/failures.txt"; continue; }
  # normalise: ISEScan's result file lands somewhere under out/ named after the
  # input; collect it to one flat place so stage 60 can glob it.
  found=$(find "$work/out" -name "*.tsv" -o -name "*.csv" | head -1)
  if [ -n "$found" ]; then cp -f "$found" "$OUT/tsv/$s.tsv"; else
    echo "NO_TABLE $s" >> "$OUT/failures.txt"; fi
done
echo "task $I done"
