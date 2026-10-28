#!/usr/bin/env bash
# Rebuild the paper's tables and figures from the archived results.
# tabulate_all_aligned.R expects the folder layout of the original project;
# input_map.tsv maps each input it reads to the archived copy in this repository.
# Usage: bash tabulation/run_tabulation.sh [output folder]   (default: current folder)
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/.." && pwd)"
OUT="${1:-$PWD}"; mkdir -p "$OUT"
RS=""
for c in "$(command -v Rscript 2>/dev/null || true)" /Library/Frameworks/R.framework/Resources/bin/Rscript \
         /opt/homebrew/bin/Rscript /usr/local/bin/Rscript /usr/bin/Rscript; do
  if [ -n "$c" ] && [ -x "$c" ]; then RS="$c"; break; fi
done
[ -n "$RS" ] || { echo "Rscript not found; install R or add it to PATH." >&2; exit 1; }
STAGE="$(mktemp -d)"; trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/codes"
TAB="$(printf '\t')"
while IFS="$TAB" read -r dest src || [ -n "$dest" ]; do
  case "$dest" in ''|\#*) continue ;; esac
  mkdir -p "$STAGE/$(dirname "$dest")"; ln -s "$ROOT/$src" "$STAGE/$dest"
done < "$HERE/input_map.tsv"
cd "$OUT" && REPRO_ROOT="$STAGE" "$RS" "$HERE/tabulate_all_aligned.R"
echo "Tables and figures written under $OUT/results/"
