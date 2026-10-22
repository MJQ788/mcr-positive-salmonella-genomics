#!/usr/bin/env bash

set -euo pipefail

PROTEINS="06_TA_systems/01_plasmid_proteins/mcr139_complete_plasmid_proteins.faa"
QUERY="06_TA_systems/02_TADB2_reference/processed/TADB2_typeII_components.faa"
DB="06_TA_systems/03_blast/database/mcr139_complete_plasmid_proteins"
RESULTS="06_TA_systems/03_blast/results/TADB2_vs_mcr139_proteins.raw.tsv"

mkdir -p 00_admin 06_TA_systems/03_blast/database 06_TA_systems/03_blast/results
test -s "$PROTEINS" || { echo "Missing $PROTEINS" >&2; exit 1; }
test -s "$QUERY" || { echo "Missing $QUERY" >&2; exit 1; }

makeblastdb \
  -in "$PROTEINS" \
  -dbtype prot \
  -parse_seqids \
  -title "937 complete mcr-positive plasmid proteins" \
  -out "$DB" \
  2>&1 | tee 00_admin/makeblastdb_mcr139_TA_proteins.log

blastp \
  -query "$QUERY" \
  -db "$DB" \
  -evalue 1e-5 \
  -seg no \
  -max_hsps 1 \
  -max_target_seqs 5000 \
  -num_threads 48 \
  -outfmt '6 qseqid sseqid pident length qlen qcovhsp qstart qend sstart send slen evalue bitscore' \
  -out "$RESULTS" \
  2>&1 | tee 00_admin/blastp_TADB2_vs_mcr139.log

