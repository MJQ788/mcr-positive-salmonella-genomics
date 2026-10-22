#!/usr/bin/env bash

set -euo pipefail

ALIGNMENT="${1:-clean.core.aln}"
OUTPUT="${2:-distances.tsv}"
snp-dists "$ALIGNMENT" > "$OUTPUT"
echo "SNP distance matrix: $OUTPUT"
