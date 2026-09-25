#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# Generic wrapper for ANY analysis step. Nothing runs on the login node.
#
# WHY THIS EXISTS: the staged compute (Arm A/B, Layer 2, ANI, downloads) was
# always submitted through sbatch, but the "quick" analysis steps -- dedup,
# the mobility join, the catalogue writer, validity checks, panel sampling,
# edge extraction -- were run directly in the login shell. Two of those were
# badly out of line: 99_mobility_join.py spawns `minimap2 -t 8` per genome
# across thousands of genomes, and the edge extraction used `sort -u -S 4G`.
# A login node is for ls, wc, head and squeue.
#
#   sbatch --export=ALL,SCRIPT=scripts/100_event_dedup.py,ARGS="--recon ... --out ..." \
#          --output=<log> slurm/analysis.sh
#
# ARGS is word-split deliberately so multi-value flags (--recon a b c) work;
# it must not contain paths with spaces.
#
# GLOBS: `set -f` disables pathname expansion during that word-split, so a
# pattern like 'shards/ani_*.tsv' reaches the script INTACT for Python to
# expand. Without it the submitting shell expands the pattern, ARGS carries 40
# separate paths, and a single-value argparse flag exits 2 with "unrecognized
# arguments" -- which is what happened to 102_rarefaction.py. Scripts that take
# globs must therefore accept either form; ours use nargs="+" plus glob.glob.
set -uo pipefail
P=/global/home/users/kh36969/fna_based_mgefinder_project
export PATH="$P/env/bin:/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
for t in python3 minimap2 samtools; do
    command -v "$t" >/dev/null || { echo "FATAL: $t not on PATH"; exit 4; }
done
: "${SCRIPT:?SCRIPT must be set}"
echo "host $(hostname)  job ${SLURM_JOB_ID:-none}  cpus ${SLURM_CPUS_PER_TASK:-?}"
echo "run  $SCRIPT ${ARGS:-}"
# shellcheck disable=SC2086
set -f                      # no pathname expansion while word-splitting ARGS
python3 "$P/$SCRIPT" ${ARGS:-}
rc=$?
echo "exit=$rc"
exit $rc
