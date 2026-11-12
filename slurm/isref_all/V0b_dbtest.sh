#!/bin/bash
# V0-4 tool self-test 2: can blastdbcmd verify a database against a manifest
# when sequence ids contain a pipe?  BLAST may read "|" as a seqid separator.
# V2-5 relies on two behaviours of makeblastdb WITHOUT -parse_seqids (the
# exact flags C_sweep.sh and V2 use):
#   * `blastdbcmd -info` reports 3 sequences and 650 residues for this input;
#   * `blastdbcmd -entry all -outfmt %t` returns the full ids ACC|CONTIG.
# Both are ASSERTED here; any deviation exits non-zero.  The -parse_seqids
# variant is run for the record only (it is not used anywhere).
#
# Usage: V0b_dbtest.sh <outdir>
# Called only by V0_prep.sh (no #SBATCH header: it is not a job of its own).
# GIT_COMMIT must be set, as for every job step (G-4); the G-4 fields and every
# tool command are written to the log, and commands are appended to
# $V0B_CMDS when V0_prep.sh passes it.
set -euo pipefail
OUT="${1:?usage: V0b_dbtest.sh <outdir>}"
: "${GIT_COMMIT:?GIT_COMMIT not set (G-4); V0b_dbtest.sh runs only inside V0_prep.sh}"
REPO=$(cd "$(dirname "$0")/../.." && pwd)
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
TD=$(mktemp -d "${TMPDIR:-/tmp}/v0b.XXXXXX")
trap 'rm -rf "$TD"' EXIT
mkdir -p "$OUT"
LOG="$OUT/v0b_dbtest.log"
: > "$LOG"
log() { echo "$*" | tee -a "$LOG"; }
CMDS="${V0B_CMDS:-$OUT/v0b_commands.txt}"
run_logged() { { printf '%q ' "$@"; echo; } | tee -a "$CMDS" >> "$LOG"; "$@"; }
# sha_of <file>: prints the sha256 or returns non-zero.  Every value that is
# written to a record is computed by a checked assignment first, never inside
# an echo, so a failing sha256sum/git cannot be masked by echo's exit status.
sha_of() { local h; h=$(sha256sum "$1" | awk '{print $1}') || return 1
           [[ "$h" =~ ^[0-9a-f]{64}$ ]] || return 1; printf '%s' "$h"; }
die() { echo "FATAL: $*"; exit 1; }
git_runtime=$(git -C "$REPO" rev-parse HEAD) || die "git rev-parse"
git_status=$(git -C "$REPO" status --porcelain --untracked-files=all) || die "git status"
if [ -n "$git_status" ]; then git_dirty=true; else git_dirty=false; fi
host=$(hostname) || die "hostname"
v_blastn=$("$BL/blastn" -version | head -1) || die "blastn -version"
v_mkdb=$("$BL/makeblastdb" -version | head -1) || die "makeblastdb -version"
v_dbcmd=$("$BL/blastdbcmd" -version | head -1) || die "blastdbcmd -version"
self_sha=$(sha_of "$0") || die "sha256 of $0"
log "git_commit_submit=$GIT_COMMIT git_commit_runtime=$git_runtime git_dirty=$git_dirty"
log "slurm_job_id=${SLURM_JOB_ID:-none} array_task=${SLURM_ARRAY_TASK_ID:-none} host=$host"
log "blastn: $v_blastn; makeblastdb: $v_mkdb; blastdbcmd: $v_dbcmd"
log "script sha256: $self_sha"

printf '>GCA_000005845.2|U00096.3\n%s\n>GCA_000005845.2|PLASMID1\n%s\n>GCF_000009999.1|CP000001.1\n%s\n' \
  "$(head -c 300 /dev/zero | tr '\0' 'A')" "$(head -c 200 /dev/zero | tr '\0' 'C')" \
  "$(head -c 150 /dev/zero | tr '\0' 'G')" > "$TD/x.fna"
printf 'GCA_000005845.2|PLASMID1\nGCA_000005845.2|U00096.3\nGCF_000009999.1|CP000001.1\n' > "$TD/want_ids"
x_sha=$(sha_of "$TD/x.fna") || die "sha256 of x.fna"
log "input: 3 contigs, 650 bp, ids with '|'; x.fna sha256 $x_sha"

# ---- asserted: no -parse_seqids (as C_sweep.sh / V2) ----
run_logged "$BL/makeblastdb" -in "$TD/x.fna" -dbtype nucl -out "$TD/xdb" -title x > "$TD/mk.log" 2>&1
info=$(run_logged "$BL/blastdbcmd" -db "$TD/xdb" -info)
n=$(printf '%s\n' "$info" | awk '/sequences;/{gsub(/,/,""); print $1}')
b=$(printf '%s\n' "$info" | awk '/sequences;/{gsub(/,/,""); print $3}')
log "no -parse_seqids: -info sequences=$n residues=$b (want 3 / 650)"
[ "$n" = 3 ] && [ "$b" = 650 ] || { log "V0B FAIL: -info counts"; exit 1; }
run_logged "$BL/blastdbcmd" -db "$TD/xdb" -entry all -outfmt '%t' | LC_ALL=C sort > "$TD/got_t"
run_logged "$BL/blastdbcmd" -db "$TD/xdb" -entry all -outfmt '%a' | LC_ALL=C sort > "$TD/got_a"
log "no -parse_seqids: %t ->"; sed 's/^/    /' "$TD/got_t" | tee -a "$LOG"
log "no -parse_seqids: %a ->"; sed 's/^/    /' "$TD/got_a" | tee -a "$LOG"
cmp -s "$TD/got_t" "$TD/want_ids" || { log "V0B FAIL: %t does not return the full ACC|CONTIG ids"; exit 1; }

# ---- record only: -parse_seqids ----
if run_logged "$BL/makeblastdb" -in "$TD/x.fna" -dbtype nucl -out "$TD/ydb" -title y -parse_seqids > "$TD/mk2.log" 2>&1; then
  log "-parse_seqids (record only): %a ->"
  run_logged "$BL/blastdbcmd" -db "$TD/ydb" -entry all -outfmt '%a' 2>&1 | sed 's/^/    /' | tee -a "$LOG"
else
  log "-parse_seqids (record only): makeblastdb failed:"; sed 's/^/    /' "$TD/mk2.log" | tee -a "$LOG"
fi
log "V0B PASS: -info counts and %t full ids as V2-5 requires"
