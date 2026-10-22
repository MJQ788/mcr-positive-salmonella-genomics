#!/usr/bin/env bash

set -euo pipefail

ALIGNMENT="${1:-clean.core.aln}"
PREFIX="${2:-iqtree_final_tree}"
IQTREE_BIN="${IQTREE_BIN:-iqtree}"

"$IQTREE_BIN" \
  -s "$ALIGNMENT" \
  -m K3Pu+ASC+R2 \
  -B 1000 \
  -nt 3 \
  --prefix "$PREFIX" \
  > "${PREFIX}.log" 2>&1
