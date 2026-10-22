#!/usr/bin/env bash
set -euo pipefail
mkdir -p pyseer_results/filtered_nondefinition_gwas

case "${1:-}" in
  full)
    pyseer \
      --phenotypes panGWAS_input/pyseer_Hybrid_phenotypes.tsv \
      --phenotype-column Hybrid \
      --pres panGWAS_input/pyseer_gene_presence_absence_min1pct.Rtab \
      --distances panGWAS_input/pyseer_coreSNP_distances.tsv \
      --max-dimensions 4 \
      --save-m pyseer_results/Hybrid_coreSNP_MDS \
      --min-af 0.01 --max-af 0.99 --cpu 40 \
      --print-samples \
      --output-patterns pyseer_results/Hybrid_gene_patterns.txt \
      > pyseer_results/panGWAS_full_pyseer.tsv \
      2> pyseer_results/panGWAS_full_pyseer.log
    ;;
  conservative)
    pyseer \
      --phenotypes panGWAS_input/pyseer_Hybrid_phenotypes.tsv \
      --phenotype-column Hybrid \
      --pres panGWAS_input/pyseer_gene_presence_absence_min1pct.filtered_conservative.Rtab \
      --distances panGWAS_input/pyseer_coreSNP_distances.tsv \
      --max-dimensions 4 \
      --save-m pyseer_results/filtered_nondefinition_gwas/Hybrid_coreSNP_MDS.filtered_conservative \
      --min-af 0.01 --max-af 0.99 --cpu 40 \
      --print-samples \
      --output-patterns pyseer_results/filtered_nondefinition_gwas/Hybrid_gene_patterns.filtered_conservative.txt \
      > pyseer_results/filtered_nondefinition_gwas/panGWAS_filtered_conservative_pyseer.tsv \
      2> pyseer_results/filtered_nondefinition_gwas/panGWAS_filtered_conservative_pyseer.log
    ;;
  strict)
    pyseer \
      --phenotypes panGWAS_input/pyseer_Hybrid_phenotypes.tsv \
      --phenotype-column Hybrid \
      --pres panGWAS_input/pyseer_gene_presence_absence_min1pct.filtered_strict.Rtab \
      --distances panGWAS_input/pyseer_coreSNP_distances.tsv \
      --max-dimensions 4 \
      --save-m pyseer_results/filtered_nondefinition_gwas/Hybrid_coreSNP_MDS.filtered_strict \
      --min-af 0.01 --max-af 0.99 --cpu 40 \
      --print-samples \
      --output-patterns pyseer_results/filtered_nondefinition_gwas/Hybrid_gene_patterns.filtered_strict.txt \
      > pyseer_results/filtered_nondefinition_gwas/panGWAS_filtered_strict_pyseer.tsv \
      2> pyseer_results/filtered_nondefinition_gwas/panGWAS_filtered_strict_pyseer.log
    ;;
  *)
    echo "Usage: bash scripts/01_pyseer_workflow.sh {full|conservative|strict}" >&2
    exit 2
    ;;
esac
