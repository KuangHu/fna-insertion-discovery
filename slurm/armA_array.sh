#!/bin/bash
#SBATCH --job-name=armA_multicopy
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-99
#SBATCH --output=%x_%A_%a.out
# Arm A is ~1 s of CPU per complete genome, so pack many genomes per array task.
#
#   sbatch --export=FNA_LIST=/path/list.txt,OUT=/path/armA slurm/armA_array.sh
set -euo pipefail
: "${FNA_LIST:?set FNA_LIST}" ; : "${OUT:?set OUT}"
REPO=/global/home/users/kh36969/fna_based_mgefinder_project
# The `fnains` core env from setup_env.sh was never built; the working runtime
# is the repo's static binaries plus claude-env (minimap2, biopython, stdlib).
export PATH="$REPO/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
mkdir -p "$OUT"

N=$SLURM_ARRAY_TASK_COUNT ; I=$SLURM_ARRAY_TASK_ID
awk -v n="$N" -v i="$I" 'NR%n==i' "$FNA_LIST" | while read -r fna; do
  [ -s "$fna" ] || continue
  s=$(basename "$fna"); s=${s%.gz}; s=${s%.fna}; s=${s%.fa}; s=${s%.fasta}
  [ -s "$OUT/${s}_families.tsv" ] && continue
  work="$fna"
  if [[ "$fna" == *.gz ]]; then
    work="$OUT/$s.fna"; zcat "$fna" > "$work"
  fi
  python3 "$REPO/scripts/20_armA_multicopy.py" "$work" \
      --out "$OUT/$s" --sample-id "$s" --threads "${SLURM_CPUS_PER_TASK:-8}" \
      >> "$OUT/$s.log" 2>&1 || echo "FAIL $s" >> "$OUT/failures.txt"
  [[ "$fna" == *.gz ]] && rm -f "$work" "$work.fai"
done
echo "task $I done"
