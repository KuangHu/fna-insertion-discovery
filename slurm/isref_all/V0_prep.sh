#!/bin/bash
#SBATCH --job-name=isall_V0
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=60G
#SBATCH --time=03:00:00
#SBATCH --exclude=n0052.cf1
#SBATCH --output=/global/scratch/users/kh36969/fna_ins_discovery/isref_all/v0/logs/%x_%j.out
# SPEC item 0 (V0 trusted manifest).  Submit from the repo root as
#   mkdir -p $W/v0/logs
#   sbatch --exclude=n0052.cf1 --export=ALL,GIT_COMMIT=$(git rev-parse HEAD) slurm/isref_all/V0_prep.sh
# n0052 is excluded (SPEC V0-7): job 27786542 failed every read there with
# Input/output error (docs/lrc_io_incident.md).  248 re-reads a genome after a
# transient EIO/ENOTCONN/ESHUTDOWN (V0-7) and fails the run at 5 such errors.
#
# Writes ONLY into a new run directory $W/v0/runs/<SLURM_JOB_ID>/ (refuses an
# existing one).  Reads, never writes, everything else: the query files, the
# corpus table and the 40,541 genome FASTA files.  It does not touch any
# legacy directory (copies/, hits/, shardlog/, verify/).
#
#   1. V0-4 tool self-test 1: outfmt 7 accounts for every query, hit or not
#      (1 hit + 2 no-hit queries; all three announced, trailer
#      "processed 3 queries").  Asserted.
#   2. V0-4 tool self-test 2: V0b_dbtest.sh (blastdbcmd -info / %t).  Asserted.
#   3. 248 build: query literals (V0-1), partition (V0-2), genome trusted
#      manifest with per-contig lengths, all recomputed from the sources
#      (V0-3), re-reading a genome after a transient I/O error (V0-7; every
#      genome's re-reads go to genome_reads.tsv), then an independent re-check
#      of what was written.  248 also
#      refuses to run unless 248, this script and V0b_dbtest.sh are
#      byte-identical to GIT_COMMIT, and records the self-test files
#      (commands, raw output, verdicts) by sha256 in v0_result.json.
# The run directory is complete only if it contains v0_result.json (written
# last by 248).  The candidate config isref_all_trusted.tsv is reviewed by the
# user before it is committed as config/isref_all_trusted.tsv (V0-5).
set -Eeuo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
V=/global/scratch/users/kh36969/fna_ins_discovery/isref_v1
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
REPO=/global/home/users/kh36969/fna_based_mgefinder_project

: "${GIT_COMMIT:?GIT_COMMIT not set: sbatch --export=ALL,GIT_COMMIT=\$(git rev-parse HEAD)}"
[ -n "${SLURM_JOB_ID:-}" ] || { echo "FATAL: run under SLURM (login node only runs small tests)"; exit 1; }
RUN="$W/v0/runs/$SLURM_JOB_ID"
mkdir -p "$W/v0/runs"
mkdir "$RUN" || { echo "FATAL: run directory exists: $RUN"; exit 1; }
# From here on, ANY non-zero exit (a version command, a self-test, V0b, 248)
# leaves $RUN/FAILED with the failing line and command, the node name and the
# SLURM job id (V0-7).  When 248 itself failed, its own record is
# $RUN/manifest/FAILED and $RUN/FAILED points to it.
LAST_LINE=?; LAST_CMD=?; TD=
trap 'LAST_LINE=$LINENO; LAST_CMD=$BASH_COMMAND' ERR
on_exit() {
  local rc=$?
  [ -n "$TD" ] && rm -rf "$TD"
  if [ "$rc" -ne 0 ] && [ ! -e "$RUN/FAILED" ]; then
    {
      echo "V0_prep.sh exit $rc at line $LAST_LINE: $LAST_CMD"
      [ -e "$RUN/manifest/FAILED" ] && echo "248 failed; its record: $RUN/manifest/FAILED"
      echo "host=$(hostname) slurm_job_id=$SLURM_JOB_ID git_commit_submit=$GIT_COMMIT"
    } > "$RUN/FAILED"
  fi
  exit "$rc"
}
trap on_exit EXIT
TT="$RUN/tooltests"; mkdir "$TT"
TD=$(mktemp -d "${TMPDIR:-/tmp}/v0.XXXXXX")
# sha_of <file>: prints the sha256 or returns non-zero.  Every value that is
# written to a record is computed by a checked assignment first, never inside
# an echo, so a failing sha256sum/git cannot be masked by echo's exit status.
sha_of() { local h; h=$(sha256sum "$1" | awk '{print $1}') || return 1
           [[ "$h" =~ ^[0-9a-f]{64}$ ]] || return 1; printf '%s' "$h"; }
die() { echo "FATAL: $*"; exit 1; }
host=$(hostname) || die "hostname"
git_runtime=$(git -C "$REPO" rev-parse HEAD) || die "git rev-parse"
echo "V0 run $SLURM_JOB_ID on $host; git_commit_submit=$GIT_COMMIT; runtime=$git_runtime"
"$BL/blastn" -version | head -1
"$BL/makeblastdb" -version | head -1
"$PY" --version

echo "############ TEST 1: outfmt 7 per-query accounting ############"
# a real 1 kb slice of the first corpus genome (must hit), a shuffled copy and
# a low-complexity sequence (must not hit)
G=$(awk -F'\t' 'NR==2{print $3}' "$V/corpus_v2.tsv")
"$PY" -I - "$G" "$TD/t_db.fna" "$TD/t_q.fna" <<'PYEOF'
import sys, random
g, dbo, qo = sys.argv[1], sys.argv[2], sys.argv[3]
seq, name = [], None
with open(g) as f:
    for line in f:
        if line.startswith(">"):
            if name: break
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
s = "".join(seq)
chunk = s[20000:21000]
assert len(chunk) == 1000, len(chunk)
with open(dbo, "w") as f:
    f.write(f">{name}\n{s[:200000]}\n")
rnd = list(chunk); random.Random(7).shuffle(rnd)
with open(qo, "w") as f:
    f.write(f">Q_hit\n{chunk}\n>Q_nohit\n{''.join(rnd)}\n>Q_nohit2\n{'ACGT'*250}\n")
PYEOF
# every tool command of the self-tests is written to commands.txt verbatim
# (G-4 full parameter strings) before it runs; 248 records the file by sha256
CMDS="$TT/commands.txt"
run_logged() { printf '%q ' "$@" >> "$CMDS"; echo >> "$CMDS"; "$@"; }
g_sha=$(sha_of "$G") || die "sha256 of test genome $G"
q_sha=$(sha_of "$TD/t_q.fna") || die "sha256 of t_q.fna"
d_sha=$(sha_of "$TD/t_db.fna") || die "sha256 of t_db.fna"
echo "# V0 job $SLURM_JOB_ID; test genome $G (sha256 $g_sha)" > "$CMDS"
echo "# TEST 1 test inputs: t_q.fna sha256 $q_sha, t_db.fna sha256 $d_sha" >> "$CMDS"
run_logged "$BL/makeblastdb" -in "$TD/t_db.fna" -dbtype nucl -out "$TD/t_db" -title t > /dev/null
run_logged "$BL/blastn" -task dc-megablast -query "$TD/t_q.fna" -db "$TD/t_db" \
  -evalue 1e-30 -perc_identity 80 -max_target_seqs 1000000 -num_threads 4 \
  -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  -out "$TT/outfmt7_test.tsv"
"$PY" -I - "$TT/outfmt7_test.tsv" <<'PYEOF' | tee "$TT/outfmt7_test.verdict"
import sys
blocks, cur, trailer, after_trailer = {}, None, None, 0
for line in open(sys.argv[1]):
    line = line.rstrip("\n")
    if trailer is not None and line and not line.startswith("#"):
        after_trailer += 1
    if line.startswith("# Query: "):
        cur = line[len("# Query: "):].split()[0]
        assert cur not in blocks, f"duplicate announcement {cur}"
        blocks[cur] = {"declared": None, "rows": 0}
    elif line.startswith("# BLAST processed "):
        trailer = line
    elif line.endswith(" hits found") and line.startswith("# "):
        blocks[cur]["declared"] = int(line.split()[1])
    elif line and not line.startswith("#"):
        q = line.split("\t")[0]
        assert q == cur, f"row for {q} under block {cur}"
        blocks[cur]["rows"] += 1
ok = (set(blocks) == {"Q_hit", "Q_nohit", "Q_nohit2"}
      and all(b["declared"] == b["rows"] for b in blocks.values())
      and blocks["Q_hit"]["rows"] >= 1
      and blocks["Q_nohit"]["rows"] == 0 and blocks["Q_nohit2"]["rows"] == 0
      and trailer == "# BLAST processed 3 queries" and after_trailer == 0)
print("blocks:", blocks)
print("trailer:", trailer)
print("V0 TEST1 PASS: all 3 queries announced, counts match, trailer processed 3" if ok
      else "V0 TEST1 FAIL")
sys.exit(0 if ok else 1)
PYEOF

echo "############ TEST 2: blastdbcmd id / count behaviour (V0b) ############"
V0B_CMDS="$CMDS" bash "$REPO/slurm/isref_all/V0b_dbtest.sh" "$TT"

echo "############ BUILD: trusted manifest (248) ############"
# SPEC V2-2 literals, the V2-4(b) partition and the G-5 BLAST version;
# 248 fails unless every one is equal
cat > "$RUN/expect.json" <<'JSON'
{"all_queries_fna_sha256": "8d5bbd5a99988d34d94f8117dc712a5e0ba23aafa40701dccea661355d1e504a",
 "query_ids_sorted_sha256": "45dd8d78c6a990be93252741f4c76592126ce3aac9605ea439d84c87ff423caa",
 "n_queries": 63095,
 "all_queries_fna_bytes": 96988793,
 "n_genomes": 40541,
 "n_shards": 64,
 "shard_sizes": "35x633,29x634",
 "blast_version": "2.17.0"}
JSON
"$PY" -I "$REPO/scripts/248_trusted_manifest.py" build \
  --queries-fna "$W/all_queries.fna" --queries-tsv "$W/all_queries.tsv" \
  --corpus "$V/corpus_v2.tsv" --expect "$RUN/expect.json" \
  --blast-bin "$BL" --repo "$REPO" --threads 64 --tooltests "$TT" \
  --code-files scripts/248_trusted_manifest.py slurm/isref_all/V0_prep.sh slurm/isref_all/V0b_dbtest.sh \
  --out "$RUN/manifest"
echo "V0 COMPLETE: $RUN"
