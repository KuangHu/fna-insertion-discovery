#!/bin/bash
# RUN THIS ON A TRANSFER / LOGIN-XFER NODE, NOT VIA sbatch.
#
# Downloads are pure network I/O. cf1 is OverSubscribe=EXCLUSIVE, so a download
# submitted there is charged a WHOLE 64-CORE NODE while using ~0 CPU. Measured:
# the K. pneumoniae download cost 225 CPU-hours and A. baumannii 40, for work
# that needs no cores at all.
#
#   ssh lrc-xfer.lbl.gov          # or your site's transfer host
#   cd /global/home/users/kh36969/fna_based_mgefinder_project
#   ./tools/download_on_transfer_node.sh saureus            # one species
#   ./tools/download_on_transfer_node.sh                    # everything queued
#
# Resumable: files already present and non-empty are skipped, so re-running
# after an interruption costs nothing. Safe to run several species in separate
# shells.
set -uo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery
Q=$W/download_queue
PAR=${PAR:-4}          # concurrent transfers. 4 is polite to NCBI; 128 drew
                       # rate-limiting and a 1.2% transient failure rate.

fetch() {
  acc="$1"; url="$2"; dest="$3/${acc}.fna"
  [ -s "$dest" ] && return 0
  tmp="$dest.part"
  for try in 1 2 3 4 5 6; do
    if curl -fsSL --max-time 600 --connect-timeout 30 -o "$tmp" "$url" \
       && gunzip -t "$tmp" 2>/dev/null; then
      gunzip -c "$tmp" > "$dest" && rm -f "$tmp" && return 0
    fi
    rm -f "$tmp"; sleep $(( (1 << try) + RANDOM % 5 ))
  done
  echo "FAIL $acc $url" >&2; return 1
}
export -f fetch

species=${1:-}
list=$( [ -n "$species" ] && echo "$Q/$species.tsv" || ls "$Q"/*.tsv )
for f in $list; do
  sp=$(basename "$f" .tsv)
  out=$W/genomes_$sp
  mkdir -p "$out"
  want=$(wc -l < "$f")
  echo "=== $sp : $want genomes -> $out"
  awk -F'\t' -v o="$out" '{print $1"\t"$2"\t"o}' "$f" \
    | xargs -P "$PAR" -n 3 bash -c 'fetch "$0" "$1" "$2"' 2> "$W/download_queue/${sp}.fail"
  have=$(ls "$out"/*.fna 2>/dev/null | wc -l)
  nf=$(grep -c '^FAIL' "$W/download_queue/${sp}.fail" 2>/dev/null || true)
  echo "    got $have / $want   failures ${nf:-0}"
  # A missing genome is INVISIBLE downstream -- no stage checks pool
  # completeness, so a partial download silently becomes a smaller species.
  if [ "$have" -ne "$want" ]; then
    echo "    INCOMPLETE -- re-run this script for $sp; it skips what is present"
  else
    echo "    COMPLETE"
  fi
done
