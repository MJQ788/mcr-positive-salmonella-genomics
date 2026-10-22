#!/usr/bin/env bash

set -euo pipefail

: "${DB_PREFIX:?Export DB_PREFIX pointing to your locally built RefSeq plasmid BLAST DB.}"
QUERY="01_mcr_reference/mcr_primary_unique.fasta"
RAW="03_mcr_screening/refseq_plasmid_primary/raw/mcr_primary_vs_refseq_plasmid.raw.tsv"

test -s "$QUERY" || { echo "Missing $QUERY" >&2; exit 1; }
mkdir -p "$(dirname "$RAW")"

blastn \
  -task blastn \
  -query "$QUERY" \
  -db "$DB_PREFIX" \
  -evalue 1e-20 \
  -perc_identity 85 \
  -dust no \
  -num_threads 48 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp qcovs qcovus stitle' \
  -out "$RAW"


