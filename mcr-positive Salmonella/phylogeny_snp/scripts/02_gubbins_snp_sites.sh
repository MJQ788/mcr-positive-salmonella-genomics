#!/usr/bin/env bash

set -euo pipefail

CORE_FULL="${1:-core.full.aln}"
CLEAN_FULL="${2:-clean.full.aln}"
GUBBINS_PREFIX="${3:-gubbins_result}"
SNP_ALIGNMENT="${4:-clean.core.aln}"
GUBBINS_THREADS="${GUBBINS_THREADS:-40}"

snippy-clean_full_aln "$CORE_FULL" > "$CLEAN_FULL"

run_gubbins.py \
  --prefix "$GUBBINS_PREFIX" \
  --threads "$GUBBINS_THREADS" \
  "$CLEAN_FULL" \
  > "${GUBBINS_PREFIX}.log" 2>&1

snp-sites -c "${GUBBINS_PREFIX}.filtered_polymorphic_sites.fasta" > "$SNP_ALIGNMENT"
echo "Recombination-filtered core-SNP alignment: $SNP_ALIGNMENT"
