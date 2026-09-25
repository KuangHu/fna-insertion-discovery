#!/bin/bash
# Mirror the Claude Code memory directory into the repo.
#
# The memory lives under ~/.claude, which is not version-controlled and has
# already been lost once (see docs/MEMORY_BACKUP.md). Run this at the end of any
# session that added or edited a memory file.
#
#   ./tools/backup_memory.sh            # memory  -> docs/memory   (backup)
#   ./tools/backup_memory.sh --restore  # docs/memory -> memory    (recover)
set -euo pipefail
REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
M="$HOME/.claude/projects/-global-home-users-kh36969-fna-based-mgefinder-project/memory"
DEST="$REPO/docs/memory"

if [ "${1:-}" = "--restore" ]; then
  [ -d "$DEST" ] || { echo "no backup at $DEST" >&2; exit 1; }
  mkdir -p "$M"
  cp -v "$DEST"/*.md "$M"/
  echo "restored $(ls -1 "$DEST"/*.md | wc -l) file(s) -> $M"
else
  [ -d "$M" ] || { echo "no memory dir at $M" >&2; exit 1; }
  mkdir -p "$DEST"
  # delete-then-copy, so a memory file removed on purpose does not linger here
  rm -f "$DEST"/*.md
  cp -v "$M"/*.md "$DEST"/
  echo "backed up $(ls -1 "$DEST"/*.md | wc -l) file(s) -> $DEST"
fi
