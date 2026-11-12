#!/bin/bash
#SBATCH --job-name=m2_corpus
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --exclude=n0052.cf1
#SBATCH --output=/global/scratch/users/kh36969/fna_ins_discovery/corpus_m2/logs/%x_%j.out
# SPEC item M2-C. Submit from the repo root, only after the M2 corpus review
# passed and the code is committed:
#   mkdir -p /global/scratch/users/kh36969/fna_ins_discovery/corpus_m2/logs
#   sbatch --export=ALL,GIT_COMMIT=$(git rev-parse HEAD) slurm/corpus/M2_corpus.sh
#
# Runs scripts/250_build_corpus.py three times into a NEW run directory
# $M/runs/<SLURM_JOB_ID>/ (an existing one is refused):
#   build/  rebuild corpus_v2.tsv from config/corpus_sources.tsv; must match
#           its pinned sha256 byte for byte (C-1)
#   gaps/   gap report against the pinned GenBank summary, cutoff 2026-09-01 (C-3)
#   pairs/  full duplicate map: reads EVERY corpus genome once (40,541 files,
#           ~180 GB; V0 read the same set in 12 min) and records each one in
#           corpus_content.tsv (C-2)
# Each sub-directory is complete only if it holds <cmd>.json (written last by
# 250). Reads, never writes, the genome files and the corpus table.
# Network: none. Every input is pinned by sha256 below.
set -Eeuo pipefail
M=/global/scratch/users/kh36969/fna_ins_discovery/corpus_m2
W=/global/scratch/users/kh36969/fna_ins_discovery
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
REPO=/global/home/users/kh36969/fna_based_mgefinder_project

CORPUS=$W/isref_v1/corpus_v2.tsv
CORPUS_SHA=b33819c772c522cc066892589ac97c7845614f49579ff250ec8b55ae585ab28e
SUMMARY=$W/ncbi_meta/assembly_summary_genbank_bacteria.20261112.txt
SUMMARY_SHA=$(cut -d' ' -f1 "$SUMMARY.sha256")
CUTOFF=2026-09-01

: "${GIT_COMMIT:?GIT_COMMIT not set: sbatch --export=ALL,GIT_COMMIT=\$(git rev-parse HEAD)}"
[ -n "${SLURM_JOB_ID:-}" ] || { echo "FATAL: run under SLURM (login node only runs small tests)"; exit 1; }
[[ "$SUMMARY_SHA" =~ ^[0-9a-f]{64}$ ]] || { echo "FATAL: no sha256 for $SUMMARY"; exit 1; }
cd "$REPO"
# the code that runs must be the code that was reviewed and committed
for f in scripts/250_build_corpus.py slurm/corpus/M2_corpus.sh config/corpus_sources.tsv; do
  git show "$GIT_COMMIT:$f" | cmp -s - "$f" \
    || { echo "FATAL: $f differs from GIT_COMMIT $GIT_COMMIT"; exit 1; }
done
RUN="$M/runs/$SLURM_JOB_ID"
mkdir -p "$M/runs"
mkdir "$RUN" || { echo "FATAL: run directory exists: $RUN"; exit 1; }
trap 'rc=$?; [ $rc -ne 0 ] && echo "M2_corpus.sh exit $rc at line $LINENO host=$(hostname) job=$SLURM_JOB_ID commit=$GIT_COMMIT" > "$RUN/FAILED"; exit $rc' EXIT
echo "M2 corpus run $SLURM_JOB_ID on $(hostname); commit $GIT_COMMIT"

"$PY" scripts/250_build_corpus.py build --sources config/corpus_sources.tsv \
  --expect-sha256 "$CORPUS_SHA" --out "$RUN/build"
"$PY" scripts/250_build_corpus.py gaps --sources config/corpus_sources.tsv \
  --corpus "$CORPUS" --corpus-sha256 "$CORPUS_SHA" \
  --summary "$SUMMARY" --summary-sha256 "$SUMMARY_SHA" --cutoff "$CUTOFF" --out "$RUN/gaps"
"$PY" scripts/250_build_corpus.py pairs \
  --corpus "$CORPUS" --corpus-sha256 "$CORPUS_SHA" \
  --summary "$SUMMARY" --summary-sha256 "$SUMMARY_SHA" \
  --workers 32 --out "$RUN/pairs"
echo "complete: $RUN"
