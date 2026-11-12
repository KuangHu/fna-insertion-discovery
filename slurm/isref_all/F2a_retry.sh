#!/bin/bash
#SBATCH --job-name=isall_F2b
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=08:00:00
#SBATCH --array=0-1
#SBATCH --output=%x_%A_%a.out
# STAGE F2a -- GenBank for the 6,034 genomes the 5,133 chosen pairs need.
# Resumable: a genome is skipped only if its .gb already exists AND parses,
# so a file truncated by an interrupted write is rebuilt rather than reused.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
S=/global/home/users/kh36969/fna_based_mgefinder_project/scripts
mkdir -p "$W/gbk" "$W/gbklog"
head -c 100 "$V/corpus_v2.tsv" >/dev/null || { echo "FATAL: scratch unreadable on $(hostname)"; exit 1; }

conv() {
  local acc=$1
  local gb="$W/gbk/$acc.gb"
  if [ -s "$gb" ] && tail -c 20 "$gb" | grep -q "//"; then return 0; fi
  local path sp
  path=$(awk -F'\t' -v a="$acc" '$1==a{print $3; exit}' "$V/corpus_v2.tsv")
  sp=$(awk -F'\t' -v a="$acc" '$1==a{print $2; exit}' "$V/corpus_v2.tsv")
  [ -n "$path" ] && [ -s "$path" ] || { echo "MISSING $acc" >> "$W/gbklog/missing.txt"; return 0; }
  if ! "$PY" "$S/234_fna_to_gbk.py" --fna "$path" --out "$gb.part" \
        --organism "$sp" > "$W/gbklog/$acc.log" 2>&1; then
    echo "FAILED $acc" >> "$W/gbklog/failed.txt"; rm -f "$gb.part"; return 0
  fi
  mv "$gb.part" "$gb"
}
export -f conv; export W V PY S

T=${SLURM_ARRAY_TASK_ID:-0}; NT=${SLURM_ARRAY_TASK_COUNT:-1}
awk -v n="$NT" -v t="$T" '(NR-1)%n==t' "$W/gbk_needed.txt" \
  | xargs -r -P 32 -I{} bash -c 'conv {}'
echo "task $T done; gbk present: $(ls "$W"/gbk/*.gb 2>/dev/null | wc -l)/$(wc -l < "$W/gbk_needed.txt")"
