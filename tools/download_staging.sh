#!/bin/bash
# Download a queue of NCBI assemblies into a STAGING directory, verified.
#
#   tools/download_staging.sh <queue.tsv> <staging_dir>
#
# queue.tsv: header line, then species<TAB>accession<TAB>...<TAB>ftp_path
# (columns found by name: species, accession, ftp_path). Each assembly lands
# at <staging_dir>/<species>/<accession>.fna.
#
# Staging is deliberately NOT one of the corpus source directories
# (scripts/250_build_corpus.py): a download never changes the corpus until a
# later, explicit decision adds the staging directory as a source.
#
# Network I/O only: run on a login or transfer node, never via sbatch (cf1 is
# EXCLUSIVE; a download there is billed a whole 64-core node). Resumable:
# an assembly whose .fna is present and whose .ok record is complete (see
# ok_valid) is skipped.
#
# Each file is verified twice before it is accepted: the .gz md5 must match
# NCBI's md5checksums.txt for that assembly (fetched with its own exit status
# checked; exactly one line must name the file), and gunzip -t must pass. Every
# step after that (sha256, size, the .ok record, the renames) is checked; the
# .ok record (md5, .fna sha256 and byte count, UTC time) is moved into place
# LAST, so a .fna without a valid .ok is never counted. Failures go to
# <staging_dir>/download.fail, one line each.
#
# Completion is judged per queued assembly, by name: each must have its .fna
# and a complete .ok. The script exits non-zero if any queued assembly
# lacks them, if download.fail is not empty, or if any download worker failed;
# files in the staging directory that are not in the queue count for nothing.
set -uo pipefail
Q=${1:?queue.tsv}
D=${2:?staging dir}
PAR=${PAR:-4}          # 4 is polite to NCBI; 128 drew rate limiting
TRIES=${TRIES:-6}      # attempts per assembly
BACKOFF=${BACKOFF:-1}  # 1 = sleep 2^try + jitter between attempts; 0 = no sleep (tests)
export TRIES BACKOFF
[[ "$PAR" =~ ^[1-9][0-9]*$ && "$TRIES" =~ ^[1-9][0-9]*$ && "$BACKOFF" =~ ^[0-9]+$ ]] \
  || { echo "FATAL: PAR and TRIES must be positive integers, BACKOFF a non-negative integer"; exit 2; }
[ -s "$Q" ] || { echo "FATAL: queue $Q is missing or empty"; exit 2; }
mkdir -p "$D" || { echo "FATAL: cannot create $D"; exit 2; }
# Parse and validate the queue ONCE, into a checked list that both the
# downloads and the final per-assembly check read. Every row must have a
# safe species label (it becomes a directory name), an NCBI accession (a file
# name) and an http(s) ftp_path with a host name (optional :port) and a path
# whose last component is the accession or starts with "<accession>_" (the
# same URL rule as scripts/250_build_corpus.py FTP); duplicates are refused.
# Any bad row, or a failure of the parser itself, stops the run before
# anything is fetched.
LIST="$D/.queue.checked"
awk -F'\t' -v listn="$LIST.n" '
  function bad(msg) { print "FATAL: queue line " NR ": " msg > "/dev/stderr"; err = 1; exit 3 }
  NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
            if (!c["species"] || !c["accession"] || !c["ftp_path"]) {
              print "FATAL: queue header must name species, accession and ftp_path" > "/dev/stderr"
              err = 1; exit 3 }
            next }
  { sp = $c["species"]; acc = $c["accession"]; ftp = $c["ftp_path"]
    if (sp !~ /^[a-z][a-z0-9_]*$/) bad("species label " sp)
    if (acc !~ /^GC[AF]_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][.][0-9]+$/) bad("accession " acc)
    if (ftp !~ /^https?:\/\/[A-Za-z0-9]([A-Za-z0-9.-]*[A-Za-z0-9])?(:[0-9]+)?\/[A-Za-z0-9._~\/:+-]*$/) bad("ftp_path " ftp)
    base = ftp; sub(/\/+$/, "", base); sub(/.*\//, "", base)
    if (base != acc && index(base, acc "_") != 1) bad("ftp_path does not name " acc)
    if (seen[acc]++) bad("duplicate accession " acc)
    print sp "\t" acc "\t" ftp; n++ }
  END { if (!err) print n + 0 > listn }' "$Q" > "$LIST"
prc=$?
want=$(cat "$LIST.n" 2>/dev/null)
nl=$(wc -l < "$LIST" 2>/dev/null)
# accepted only when every value was measured and agrees (a failed
# measurement can never pass for agreement)
if [ "$prc" = 0 ] && [[ "$want" =~ ^[0-9]+$ && "$nl" =~ ^[0-9]+$ ]] && [ "$nl" -eq "$want" ]; then
  :
else
  echo "FATAL: queue $Q rejected (parser exit $prc); nothing downloaded"
  rm -f "$LIST" "$LIST.n"
  exit 2
fi
: > "$D/download.fail" || { echo "FATAL: cannot write $D/download.fail"; exit 2; }

ok_valid() {      # ok_valid <ok> <fna> <accession>: the complete record of this .fna
  # exactly the six lines fetch() writes, each once: accession (this one),
  # url, gz_md5 (32 hex), fna_sha256 (64 hex), fna_bytes (= the .fna's size
  # now), utc (YYYY-MM-DDTHH:MM:SSZ)
  local size
  [ -s "$1" ] && [ -s "$2" ] || return 1
  size=$(stat -c %s "$2") && [[ "$size" =~ ^[1-9][0-9]*$ ]] || return 1
  awk -F'\t' -v acc="$3" -v size="$size" '
    function hex(v, n) { return length(v) == n && v !~ /[^0-9a-f]/ }
    NF != 2 || seen[$1]++ { bad = 1; next }
    $1 == "accession"  { bad = bad || $2 != acc; next }
    $1 == "url"        { bad = bad || $2 !~ /^https?:\/\/[^ ]+_genomic\.fna\.gz$/; next }
    $1 == "gz_md5"     { bad = bad || !hex($2, 32); next }
    $1 == "fna_sha256" { bad = bad || !hex($2, 64); next }
    $1 == "fna_bytes"  { bad = bad || $2 !~ /^[1-9][0-9]*$/ || $2 != size; next }
    $1 == "utc"        { bad = bad || $2 !~ /^[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]Z$/; next }
                       { bad = 1 }
    END { exit (bad || NR != 6) ? 1 : 0 }' "$1"
}

fetch() {
  sp="$1"; acc="$2"; ftp="${3%/}"; D="$4"
  out="$D/$sp"; dest="$out/$acc.fna"; ok="$dest.ok"
  ok_valid "$ok" "$dest" "$acc" && return 0
  base=$(basename "$ftp")
  url="$ftp/${base}_genomic.fna.gz"
  tmp="$dest.gz.part"; sums="$dest.md5.part"
  if mkdir -p "$out"; then
    for try in $(seq 1 "$TRIES"); do
      # nothing of a failed or unfinished earlier attempt survives
      rm -f "$tmp" "$sums" "$dest.part" "$ok.part" "$dest" "$ok"
      # NCBI's md5 list goes to a file so that curl's own exit status counts;
      # exactly one line must name the file
      if curl -fsSL --max-time 120 --connect-timeout 30 -o "$sums" "$ftp/md5checksums.txt" \
         && want=$(awk -v f="${base}_genomic.fna.gz" '$2 == f || $2 == "./" f {print $1}' "$sums") \
         && [[ "$want" =~ ^[0-9a-f]{32}$ ]] \
         && curl -fsSL --max-time 600 --connect-timeout 30 -o "$tmp" "$url" \
         && got=$(md5sum < "$tmp") && [ "${got%% *}" = "$want" ] \
         && gunzip -t "$tmp" 2>/dev/null \
         && gunzip -c "$tmp" > "$dest.part" \
         && fsha=$(sha256sum < "$dest.part") && fsha=${fsha%% *} \
         && [[ "$fsha" =~ ^[0-9a-f]{64}$ ]] \
         && fbytes=$(stat -c %s "$dest.part") && [[ "$fbytes" =~ ^[1-9][0-9]*$ ]] \
         && now=$(date -u +%Y-%m-%dT%H:%M:%SZ) \
         && printf 'accession\t%s\nurl\t%s\ngz_md5\t%s\nfna_sha256\t%s\nfna_bytes\t%s\nutc\t%s\n' \
              "$acc" "$url" "$want" "$fsha" "$fbytes" "$now" > "$ok.part" \
         && rm -f "$ok" && mv "$dest.part" "$dest" && mv "$ok.part" "$ok"; then
        rm -f "$tmp" "$sums"
        return 0
      fi
      rm -f "$tmp" "$sums" "$dest.part" "$ok.part" "$dest" "$ok"
      [ "$try" -lt "$TRIES" ] && sleep $(( BACKOFF * ((1 << try) + RANDOM % 5) ))
    done
  fi
  echo "FAIL $sp $acc $url" >> "$D/download.fail"
  return 1
}
export -f fetch ok_valid

export D
xrc=0
if [ "$want" -gt 0 ]; then
  xargs -P "$PAR" -d '\n' -n 1 bash -c \
    'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; fetch "$a" "$b" "$c" "$D"' < "$LIST"
  xrc=$?
fi

seen=0; missing=0
while IFS=$'\t' read -r sp acc ftp; do
  seen=$((seen + 1))
  if ! ok_valid "$D/$sp/$acc.fna.ok" "$D/$sp/$acc.fna" "$acc"; then
    missing=$((missing + 1))
    echo "MISSING $sp $acc"
  fi
done < "$LIST"
nf=$(wc -l < "$D/download.fail") && [[ "$nf" =~ ^[0-9]+$ ]] \
  || { echo "FATAL: cannot read $D/download.fail; completion not established"; exit 1; }
echo "queue $want, verified $((want - missing)), failures $nf, worker exit $xrc"
# COMPLETE only when every count was measured and is what completion needs;
# any other outcome, including a failed test, is INCOMPLETE
if [ "$seen" -eq "$want" ] && [ "$missing" -eq 0 ] && [ "$nf" -eq 0 ] && [ "$xrc" -eq 0 ]; then
  echo COMPLETE
  exit 0
fi
[ "$seen" -eq "$want" ] || echo "FATAL: checked $seen of $want queued assemblies"
echo "INCOMPLETE: re-run; verified assemblies are skipped"
exit 1
