#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="${INPUT_DIR:-genomes}"
OUTPUT_DIR="${OUTPUT_DIR:-snippy_output}"
REFERENCE="${REFERENCE:-reference/re.gbk}"
SNIPPY_CPUS="${SNIPPY_CPUS:-50}"
CORE_PREFIX="${CORE_PREFIX:-core}"

[[ -d "$INPUT_DIR" ]] || { echo "Input genome directory missing: $INPUT_DIR" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference file missing: $REFERENCE" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR"
shopt -s nullglob
genomes=( "$INPUT_DIR"/*.fna )
(( ${#genomes[@]} > 0 )) || { echo "No .fna assemblies found in $INPUT_DIR" >&2; exit 1; }

sample_dirs=()
for genome in "${genomes[@]}"; do
    sample_name="$(basename "$genome" .fna)"
    echo "Processing: $sample_name"
    outdir="$OUTPUT_DIR/$sample_name"
    snippy --cpus "$SNIPPY_CPUS" \
           --outdir "$outdir" \
           --ref "$REFERENCE" \
           --ctgs "$genome"
    sample_dirs+=( "$outdir" )
done

# Combine all isolate Snippy folders into core alignment.
snippy-core --prefix "$CORE_PREFIX" --ref "$REFERENCE" "${sample_dirs[@]}"
echo "Core whole-reference alignment: ${CORE_PREFIX}.full.aln"
