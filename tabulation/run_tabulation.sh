#!/usr/bin/env bash
# Rebuild the paper's tables and figures from the archived results.
# tabulate_all_aligned.R expects the folder layout of the original project;
# input_map.tsv maps each input it reads to the archived copy in this repository.
# Usage: bash tabulation/run_tabulation.sh [output folder]   (default: current folder)
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/.." && pwd)"
OUT="${1:-$PWD}"; mkdir -p "$OUT"
STAGE="$(mktemp -d)"; trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/codes"
TAB="$(printf '\t')"
while IFS="$TAB" read -r dest src || [ -n "$dest" ]; do
  case "$dest" in ''|\#*) continue ;; esac
  mkdir -p "$STAGE/$(dirname "$dest")"; ln -s "$ROOT/$src" "$STAGE/$dest"
done < "$HERE/input_map.tsv"
cd "$OUT" && REPRO_ROOT="$STAGE" Rscript "$HERE/tabulate_all_aligned.R"
echo "Tables and figures written under $OUT/results/"
