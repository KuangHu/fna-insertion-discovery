#!/bin/bash
# SPIRE MAG download. RUN BY HAND ON THE TRANSFER NODE, never through sbatch.
#
#   ssh lrc-xfer.lbl.gov
#   bash /global/home/users/kh36969/fna_based_mgefinder_project/tools/download_spire_on_transfer_node.sh
#
# WHY NOT SBATCH: cf1 is OverSubscribe=EXCLUSIVE, so a job is billed a whole
# 64-core node whatever it requests. Downloads are pure network I/O and used
# ~0 CPU -- the K. pneumoniae genome download cost 225 CPU-hours for nothing.
# The transfer node costs 0 SU.
#
# Concurrency 4 with backoff. 128 parallel connections drew rate-limiting and
# a 1.2% transient failure rate on NCBI; there is no reason to assume EMBL is
# more tolerant, and a partial pool is invisible downstream.
set -uo pipefail
S=/global/scratch/users/kh36969/spire
D=$S/dl5
OUT=$S/genomes
PAR=4
mkdir -p "$OUT"

get_one(){
  id=$1; out=$2
  f="$out/$id.fna"
  # already have it, and it is not an error page or a truncated write
  [ -s "$f" ] && head -c 1 "$f" | grep -q '>' && return 0
  for attempt in 1 2 3 4 5; do
    code=$(curl -sS -L --max-time 180 --retry 0 -w '%{http_code}' \
             -o "$f.gz" "https://spire.embl.de/download_file/$id" 2>/dev/null)
    if [ "$code" = "200" ] && gzip -t "$f.gz" 2>/dev/null; then
      gunzip -f "$f.gz" 2>/dev/null && [ -s "$f" ] && \
        head -c 1 "$f" | grep -q '>' && return 0
    fi
    rm -f "$f.gz" "$f"
    # 404 is PERMANENT -- the MAG is in the metadata but not in the file
    # store. Measured on Akkermansia: 34 of 1031 (3.3%), all consistent 404s
    # on clean single requests. Retrying one costs 110s of backoff for
    # nothing, which across the full set is hours of waiting.
    if [ "$code" = "404" ]; then
      echo "$id" >> "$out/_absent_404.txt"
      return 1
    fi
    sleep $((attempt * attempt * 2))
  done
  echo "$id" >> "$out/_failed.txt"
  return 1
}
export -f get_one

for ids in "$D"/*.ids; do
  sp=$(basename "$ids" .ids)
  o="$OUT/$sp"; mkdir -p "$o"
  want=$(wc -l < "$ids")
  echo "### $sp : $want MAGs -> $o"
  xargs -P $PAR -I{} bash -c 'get_one "$@"' _ {} "$o" < "$ids"
  got=$(ls "$o"/*.fna 2>/dev/null | wc -l)
  fail=$([ -f "$o/_failed.txt" ] && sort -u "$o/_failed.txt" | wc -l || echo 0)
  # A partial pool is invisible to every later stage -- no stage checks pool
  # completeness -- so say so here, loudly, rather than letting it pass.
  if [ "$got" -eq "$want" ]; then
    echo "    COMPLETE  got $got / want $want"
  else
    echo "    INCOMPLETE  got $got / want $want  (failed $fail)"
  fi
done

echo
echo "### totals"
for d in "$OUT"/*/; do
  printf "  %-34s %6s fna  %8s\n" "$(basename "$d")" \
    "$(ls "$d"/*.fna 2>/dev/null | wc -l)" "$(du -sh "$d" 2>/dev/null | cut -f1)"
done
