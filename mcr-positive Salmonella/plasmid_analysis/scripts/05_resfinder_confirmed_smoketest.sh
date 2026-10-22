#!/usr/bin/env bash

set -euo pipefail
ARG_DIR="04_backbone_ARG/abricate_resfinder"
DB_ROOT="$ARG_DIR/db"
DB_NAME="resfinder_eecf0aa20759"
SMOKE="$ARG_DIR/smoke_test/three_plasmids.fna"
mkdir -p 00_admin "$ARG_DIR/smoke_test"
test -s "$DB_ROOT/$DB_NAME/sequences" || {
  echo "Missing frozen database sequence file" >&2; exit 1;
}
test -s "$SMOKE" || { echo "Missing three-plasmid smoke-test input" >&2; exit 1; }

makeblastdb \
  -in "$DB_ROOT/$DB_NAME/sequences" \
  -dbtype nucl \
  -title "$DB_NAME" \
  -hash_index \
  2>&1 | tee 00_admin/abricate_resfinder_makeblastdb.log

abricate \
  --datadir "$DB_ROOT" \
  --db "$DB_NAME" \
  --minid 80 \
  --mincov 80 \
  --threads 8 \
  "$SMOKE" \
  > "$ARG_DIR/smoke_test/abricate_resfinder.tsv" \
  2> "$ARG_DIR/smoke_test/abricate_resfinder.log"
