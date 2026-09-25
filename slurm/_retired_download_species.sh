#!/bin/bash
#SBATCH --job-name=dl_species
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --account=pc_rubinlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --array=0-15
#SBATCH --output=/global/scratch/users/kh36969/fna_ins_discovery/census/logs/dl_%A_%a.out
# Downloads one species' complete+chromosome assemblies to scratch.
# NOTE: compute-node curl is 7.61.1 -- it predates --retry-all-errors (needs
# 7.71), which silently failed every file on an earlier submission. Retry is
# therefore done in the loop, not by curl.
set -u
MANIFEST=${MANIFEST:?}
OUTDIR=${OUTDIR:?}
NTASK=${SLURM_ARRAY_TASK_COUNT:-16}
TID=${SLURM_ARRAY_TASK_ID:-0}
mkdir -p "$OUTDIR"
FAILLOG="${FAILDIR:-$OUTDIR/..}/dlfail_${TID}.log"
mkdir -p "$(dirname "$FAILLOG")"
# Default concurrency 2 per task. At 16 tasks that is 32 connections; the
# previous 8 (=128 connections) drew rate-limiting from NCBI.

fetch() {
  acc="$1"; url="$2"; dest="$OUTDIR/${acc}.fna"
  [ -s "$dest" ] && return 0
  tmp="$dest.gz.part"
  for try in 1 2 3 4 5 6; do
    if curl -fsSL --max-time 600 --connect-timeout 30 -o "$tmp" "$url"; then
      if gunzip -t "$tmp" 2>/dev/null; then
        gunzip -c "$tmp" > "$dest" && rm -f "$tmp" && return 0
      fi
    fi
    rm -f "$tmp"
    # EXPONENTIAL backoff with jitter, not linear. Measured on A. baumannii:
    # 16 tasks x 8 curls = 128 concurrent connections to NCBI produced a 1.2%
    # transient failure rate (24/1956) that exhausted a linear 5-try loop.
    # All 24 URLs returned 200 on a later serial probe. The failures were
    # LOUD (FAIL lines), so they were recoverable -- silently, the species
    # would have run at 98.8% coverage with nothing downstream showing it.
    sleep $(( (1 << try) + RANDOM % 5 ))
  done
  echo "FAIL $acc $url" >&2; return 1
}
export -f fetch; export OUTDIR
awk -v t="$TID" -v n="$NTASK" 'NR % n == t' "$MANIFEST" \
  | awk -F'\t' '{print $1"\t"$3}' \
  | xargs -P "${DL_PAR:-2}" -n 2 bash -c 'fetch "$0" "$1"' 2> "$FAILLOG"
cat "$FAILLOG" >&2
nf=$(grep -c "^FAIL" "$FAILLOG" 2>/dev/null || true)
nf=${nf:-0}
echo "[task $TID] done: $(ls "$OUTDIR" | wc -l) files present, $nf failures"
# ASSERTION: a missing genome is INVISIBLE downstream -- no stage checks that the
# pool is complete, so a partial download silently becomes a smaller species.
# Measured: 24/1956 (1.2%) transient failures at 128 concurrent connections.
if [ "$nf" -gt 0 ]; then
    echo "FATAL: task $TID had $nf download failures; rerun (fetch skips files"
    echo "       already present) or lower DL_PAR before continuing"
    exit 6
fi
