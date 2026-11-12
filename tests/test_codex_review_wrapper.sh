#!/bin/bash
# Hermetic tests for tools/codex-review (v4).
#
# Uses a throwaway git repo, a throwaway $HOME with fake codex profiles and a
# stub `codex` first on PATH.  Never touches the real repository, the real
# codex binary or any real review.  Usage:
#   TMPDIR=<scratch> bash tests/test_codex_review_wrapper.sh [path/to/codex-review]
set -u
WR=$(readlink -f "${1:-$(dirname "$0")/../tools/codex-review}")
[ -x "$WR" ] || { echo "wrapper not executable: $WR"; exit 2; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/crtest.XXXXXX")
UNRELATED=""
cleanup() {
  [ -n "$UNRELATED" ] && kill "$UNRELATED" 2>/dev/null
  for f in "$TMP"/stub/*.pid; do [ -f "$f" ] && kill -KILL "$(cat "$f")" 2>/dev/null; done
  rm -rf "$TMP"
}
trap cleanup EXIT
PASS=0; FAIL=0
ok()  { echo "PASS  $1"; PASS=$((PASS + 1)); }
bad() { echo "FAIL  $1 :: ${2:-}"; FAIL=$((FAIL + 1)); }
alive() { [ -n "${1:-}" ] && kill -0 "$1" 2>/dev/null; }

# ---- fake HOME with profiles, stub codex ----
export HOME="$TMP/home"; mkdir -p "$HOME/.codex" "$TMP/stub"
printf 'model = "gpt-6.1-sol"\nmodel_reasoning_effort = "high"\n' > "$HOME/.codex/review.config.toml"
export STUB_DIR="$TMP/stub"
cat > "$TMP/stub/codex" <<'EOF'
#!/bin/bash
if [ "${1:-}" = "--version" ]; then echo "codex-cli STUB"; exit 0; fi
out=""; while [ $# -gt 0 ]; do [ "$1" = "-o" ] && { out="$2"; shift; }; shift; done
cat > /dev/null
echo $$ > "$STUB_DIR/stub.pid"
printf 'OpenAI Codex STUB\n--------\nmodel: %s\nreasoning effort: high\nreasoning summaries: %s\n--------\n' "${STUB_MODEL:-gpt-6.1-sol}" "${STUB_SUMMARY:-none}"
case "${STUB_MODE:-pass}" in
  pass)      printf 'VERDICT: PASS\n无问题。\n' > "$out" ;;
  fail)      printf 'VERDICT: FAIL\n1. x\n' > "$out" ;;
  rerror)    printf 'VERDICT: ERROR\n材料截断\n' > "$out" ;;
  silent)    : ;;
  truncated) printf 'VERDICT: PA' > "$out" ;;
  crash_after_write) printf 'VERDICT: PASS\nhalf' > "$out"; exit 3 ;;
  hang)      sleep 300 & echo $! > "$STUB_DIR/grandchild.pid"
             printf 'VERDICT: PASS\npartial-before-hang\n' > "$out"; wait ;;
  hang_ignore_term)
             trap '' TERM; sleep 300 & echo $! > "$STUB_DIR/grandchild.pid"
             printf 'VERDICT: PASS\npartial-ignoring-term\n' > "$out"
             while :; do sleep 1; done ;;
esac
exit 0
EOF
chmod +x "$TMP/stub/codex"
export PATH="$TMP/stub:$PATH"

# ---- throwaway repo ----
R="$TMP/repo"; mkdir -p "$R/docs"; cd "$R" || exit 2
git init -q; git config user.email t@t; git config user.name t
echo "# agents" > AGENTS.md; printf '# spec\nline kept\nline changed NEW\n' > docs/SPEC.md; echo "committed dep" > a.txt
printf '# spec\nline kept\nline changed OLD\nOLD-ONLY-MARKER\n' > old.txt
printf '[T item] t-item\na.txt\nb.txt\n' > docs/review_deps.txt
echo "secret.log" > .gitignore
git add -A; git commit -qm base
echo "untracked dep" > b.txt
A=docs/reviews/t-item
meta() { sed -n "s/^$2: //p" "$A/r$1.md"; }      # meta <NN> <field>
run()  { "$WR" "T item" review 2>&1; }
first_line() { printf '%s\n' "$1" | head -1; }
reset_stub() { rm -f "$STUB_DIR"/*.pid; }

echo "== material checks (fail before codex) =="
expect_err() { local name="$1" frag="$2"; shift 2; local o; o=$("$@" 2>&1)
  if [ "$(first_line "$o")" = "VERDICT: ERROR" ] && printf '%s' "$o" | grep -qF "$frag"; then ok "$name"; else bad "$name" "$(printf '%s' "$o" | head -2 | tr '\n' ' ')"; fi; }
o=$(CODEX_REVIEW_DRYRUN=1 "$WR" "T item" review 2>&1)
printf '%s' "$o" | grep -q 'DRY RUN' && grep -q '^committed dep$' .review/input.txt && ok "dry run sends committed dep in full" || bad "dry run" "$o"
[ ! -d "$A" ] && ok "dry run does not archive" || bad "dry run archived"
mv b.txt b.bak; expect_err "missing dep"   "必需材料缺失：b.txt" "$WR" "T item" review; mv b.bak b.txt
printf '\x00\x01' > b.txt; expect_err "binary dep" "二进制文件：b.txt" "$WR" "T item" review; echo "untracked dep" > b.txt
head -c 300 /dev/zero | tr '\0' x > b.txt; expect_err "oversized dep" "会被截断：b.txt" env CODEX_REVIEW_MAXB=200 "$WR" "T item" review; echo "untracked dep" > b.txt
echo x > secret.log; printf '[T item] t-item\nsecret.log\n' > docs/review_deps.txt
expect_err "gitignored dep" "被 .gitignore 忽略" "$WR" "T item" review
printf '[T item] t-item\na.txt\nb.txt\n' > docs/review_deps.txt
expect_err "undeclared item"     "出现 0 次"              "$WR" "No such" review
expect_err "model pin (profile)" "要求 other-model"       env CODEX_REVIEW_MODEL=other-model "$WR" "T item" review
expect_err "bad timeout value"   "CODEX_REVIEW_TIMEOUT"   env CODEX_REVIEW_TIMEOUT=abc "$WR" "T item" review

echo "== DIFF material (v4) =="
printf '[T item] t-item\na.txt\nb.txt\nDIFF old.txt docs/SPEC.md\n' > docs/review_deps.txt
o=$(CODEX_REVIEW_DRYRUN=1 "$WR" "T item" review 2>&1)
grep -qF '=== 差异（diff -u）: old.txt -> docs/SPEC.md ===' .review/input.txt && grep -qx -- '-line changed OLD' .review/input.txt && grep -qx '+line changed NEW' .review/input.txt \
  && ok "DIFF sends diff -u of old -> new" || bad "DIFF section" "$o"
! grep -qF '=== 必需材料: old.txt' .review/input.txt && [ "$(grep -c 'OLD-ONLY-MARKER' .review/input.txt)" = 1 ] \
  && ok "DIFF old file not sent in full (only its diff lines)" || bad "old sent in full"
grep -qF '=== 差异材料清单' .review/input.txt && grep -q '^old.txt -> docs/SPEC.md.*old tracked-clean' .review/input.txt \
  && ok "DIFF listing records both sides' state and sha256" || bad "DIFF listing"
! printf '%s\n' "$(sed -n '/=== 其他未跟踪文件/,$p' .review/input.txt)" | grep -qx 'old.txt' && ok "DIFF old file not listed as 'other untracked'" || bad "other list"
mv old.txt old.bak; expect_err "DIFF old missing" "必需材料缺失：old.txt" "$WR" "T item" review; mv old.bak old.txt
printf '[T item] t-item\na.txt\nDIFF old.txt\n' > docs/review_deps.txt
expect_err "DIFF malformed line" "DIFF 行" "$WR" "T item" review
echo x > secret.log; printf '[T item] t-item\na.txt\nDIFF secret.log docs/SPEC.md\n' > docs/review_deps.txt
expect_err "DIFF gitignored side" "被 .gitignore 忽略" "$WR" "T item" review
printf '[T item] t-item\na.txt\nb.txt\n' > docs/review_deps.txt

echo "== normal runs and round numbering =="
rm -rf "$A"; reset_stub
o=$(STUB_MODE=pass run); n=01
[ "$(first_line "$o")" = "VERDICT: PASS" ] && [ "$(meta $n review_round)" = 1 ] && [ "$(meta $n completed_round)" = yes ] \
  && ok "PASS -> review_round 1" || bad "pass numbering" "$(meta $n review_round)"
o=$(STUB_MODE=silent run); n=02
[ "$(meta $n verdict)" = ERROR ] && [ "$(meta $n review_round)" = none ] && ok "ERROR consumes no round" || bad "error round" "$(meta $n review_round)"
o=$(STUB_MODE=fail run); n=03
[ "$(meta $n review_round)" = 2 ] && ok "FAIL after an ERROR is review_round 2 (run r03)" || bad "fail numbering" "$(meta $n review_round)"
id1=$(meta 01 attempt_id); id3=$(meta 03 attempt_id)
[ -n "$id1" ] && [ "$id1" != "$id3" ] && ok "attempt ids unique ($id3)" || bad "attempt ids" "$id1 / $id3"
o=$(STUB_MODE=rerror run); n=04
[ "$(meta $n verdict)" = ERROR ] && [ "$(meta $n completed_round)" = no ] && ! grep -q '已作废' "$A/r$n.md" \
  && ok "reviewer-issued ERROR archived, not marked voided, no round" || bad "reviewer error"
o=$(STUB_MODE=pass STUB_MODEL=gpt-6-astra run); n=05
[ "$(meta $n verdict)" = ERROR ] && grep -q '实际使用的模型是「gpt-6-astra」' "$A/r$n.md" && grep -q '已作废' "$A/r$n.md" \
  && ok "observed-model mismatch voids a PASS" || bad "observed model"

echo "== stale verdict rejection =="
printf 'VERDICT: PASS\nstale\n' > .review/review.md; printf 'VERDICT: PASS\nstale\n' > .review/review_raw.md
o=$(STUB_MODE=silent run); n=06
[ "$(first_line "$o")" = "VERDICT: ERROR" ] && head -1 .review/review.md | grep -qx 'VERDICT: ERROR' && ! grep -q stale .review/review.md \
  && ok "stale PASS in .review/ is never reused (silent codex -> ERROR)" || bad "stale verdict" "$(head -2 .review/review.md | tr '\n' ' ')"
[ ! -e .review/review_raw.md ] && ok "stale review_raw.md removed at start" || bad "stale raw survives"

echo "== interrupted output =="
o=$(STUB_MODE=truncated run); n=07
printf '%s' "$o" | grep -q '首行格式不符' && [ "$(meta $n verdict)" = ERROR ] && ok "truncated 'VERDICT: PA' -> ERROR" || bad "truncated" "$o"
o=$(STUB_MODE=crash_after_write run); n=08
[ "$(meta $n verdict)" = ERROR ] && grep -q 'codex 退出码 3' "$A/r$n.md" && grep -q '^half' "$A/r$n.md" \
  && ok "PASS written then crash -> ERROR, partial output kept as voided" || bad "crash after write"

echo "== internal timeout: only this attempt's process group dies =="
sleep 300 & UNRELATED=$!
reset_stub; t0=$(date +%s)
o=$(STUB_MODE=hang CODEX_REVIEW_TIMEOUT=3 CODEX_REVIEW_KILL_AFTER=2 run); n=09; t1=$(date +%s)
sp=$(cat "$STUB_DIR/stub.pid" 2>/dev/null); gp=$(cat "$STUB_DIR/grandchild.pid" 2>/dev/null)
[ "$(first_line "$o")" = "VERDICT: ERROR" ] && grep -q '超过内部超时 3s' "$A/r$n.md" && ok "timeout -> ERROR in $((t1 - t0))s" || bad "timeout verdict" "$o"
grep -q 'partial-before-hang' "$A/r$n.md" && grep -q '已作废' "$A/r$n.md" && ok "partial output before hang preserved, voided" || bad "partial on timeout"
[ "$(meta $n verdict)" = ERROR ] && ! grep -q '^VERDICT: PASS$' .review/review.md && ok "the partial 'VERDICT: PASS' did not become the verdict" || bad "partial became verdict"
! alive "$sp" && ! alive "$gp" && ok "stub and its grandchild terminated" || bad "children survived" "stub=$(alive "$sp" && echo alive) gc=$(alive "$gp" && echo alive)"
alive "$UNRELATED" && ok "unrelated process untouched" || bad "unrelated process killed"

reset_stub
o=$(STUB_MODE=hang_ignore_term CODEX_REVIEW_TIMEOUT=2 CODEX_REVIEW_KILL_AFTER=2 run); n=10
sp=$(cat "$STUB_DIR/stub.pid" 2>/dev/null); gp=$(cat "$STUB_DIR/grandchild.pid" 2>/dev/null)
[ "$(meta $n verdict)" = ERROR ] && ! alive "$sp" && ! alive "$gp" && ok "TERM-ignoring attempt killed via KILL (codex_exit $(meta $n codex_exit))" || bad "kill-after" "stub=$(alive "$sp" && echo alive) gc=$(alive "$gp" && echo alive)"
alive "$UNRELATED" && ok "unrelated process still untouched" || bad "unrelated killed (2)"

echo "== outer TERM (the 900 s limit) =="
reset_stub
o=$(STUB_MODE=hang CODEX_REVIEW_TIMEOUT=60 CODEX_REVIEW_KILL_AFTER=2 timeout 4 "$WR" "T item" review 2>&1); n=11
sleep 1
sp=$(cat "$STUB_DIR/stub.pid" 2>/dev/null); gp=$(cat "$STUB_DIR/grandchild.pid" 2>/dev/null)
[ "$(meta $n verdict)" = ERROR ] && grep -q '收到信号 TERM' "$A/r$n.md" && ok "outer TERM -> wrapper archives ERROR" || bad "outer TERM" "$(meta $n verdict) $(meta $n error_reason)"
! alive "$sp" && ! alive "$gp" && ok "outer TERM: attempt group terminated" || bad "outer TERM children"
alive "$UNRELATED" && ok "outer TERM: unrelated process untouched" || bad "unrelated killed (3)"

echo "== outer SIGKILL: pre-launch record survives =="
reset_stub
STUB_MODE=hang CODEX_REVIEW_TIMEOUT=5 CODEX_REVIEW_KILL_AFTER=2 timeout -s KILL 2 "$WR" "T item" review >/dev/null 2>&1; n=12
[ "$(meta $n verdict)" = INCOMPLETE ] && [ "$(meta $n completed_round)" = no ] && ok "wrapper SIGKILLed: attempt left as INCOMPLETE record" || bad "incomplete record" "$(meta $n verdict)"
[ ! -e .review/review.md ] && ok "no verdict file left behind after SIGKILL" || bad "verdict file after SIGKILL" "$(head -1 .review/review.md)"
sleep 9
sp=$(cat "$STUB_DIR/stub.pid" 2>/dev/null); gp=$(cat "$STUB_DIR/grandchild.pid" 2>/dev/null)
! alive "$sp" && ! alive "$gp" && ok "orphaned attempt still bounded by its own timeout" || bad "orphan survived" "stub=$(alive "$sp" && echo alive) gc=$(alive "$gp" && echo alive)"

echo "== earlier archives are never rewritten =="
before=$(sha256sum "$A"/r0[1-9].md "$A"/r1[0-2].md | sha256sum)
o=$(STUB_MODE=pass run)
after=$(sha256sum "$A"/r0[1-9].md "$A"/r1[0-2].md | sha256sum)
[ "$before" = "$after" ] && ok "r01-r12 byte-identical after a new run (r13)" || bad "old archives changed"
[ "$(meta 13 review_round)" = 3 ] && ok "review_round continues from completed runs only (r13 -> round 3)" || bad "round after many errors" "$(meta 13 review_round)"

echo "== v4 archive records DIFF files and reasoning summary =="
printf '[T item] t-item\na.txt\nb.txt\nDIFF old.txt docs/SPEC.md\n' > docs/review_deps.txt
printf 'model = "gpt-6.1-sol"\nmodel_reasoning_effort = "high"\nmodel_reasoning_summary = "concise"\n' > "$HOME/.codex/review.config.toml"
o=$(STUB_MODE=pass STUB_SUMMARY=concise run); n=14
grep -q '^  - old.txt -> docs/SPEC.md | old tracked-clean' "$A/r$n.md" && ok "archive lists diff_files" || bad "diff_files" "$(sed -n '/diff_files/,+1p' "$A/r$n.md")"
[ "$(meta $n configured_summary)" = concise ] && [ "$(meta $n observed_summary)" = concise ] && ok "archive records configured/observed reasoning summary" || bad "summary meta" "$(meta $n configured_summary)/$(meta $n observed_summary)"
[ "$(meta $n timeout_s)" = 2400 ] && ok "default internal timeout is 2400 s" || bad "default timeout" "$(meta $n timeout_s)"

echo; echo "SUMMARY: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
