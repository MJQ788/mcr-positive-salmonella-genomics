Geographic standardization and BioProject sensitivity analyses

Contents

  R_regional_sensitivity.R
      R script for geographic standardization and BioProject sensitivity
      analyses of Hybrid proportions.

  TableS1_for_R.csv
      Isolate-level metadata table with 1,514 rows.

Requirements

  R 4.5.2
  Base R only; no additional R packages required.

Input

  TableS1_for_R.csv must contain the following columns:

  BioProject
      BioProject identifier.

  Continent
      Continent label.

  Country
      Country label.

  MLST
      MLST sequence type.

  mcr subtype
      mcr subtype label.

  Isolation source
      Source label; allowed values are Human, Chicken, Environment, and Egg.

  Risk_group
      Risk group label; Hybrid is coded as 1, otherwise 0.

  Group is assigned as Human when Isolation source is Human, and as
  Chicken_related otherwise.

  The script includes fixed checks for the supplied dataset:
  1,514 rows; 879 Human; 635 Chicken_related; 223 Hybrid;
  1,442 rows in Asia, North America, and Europe; 202 Hybrid in these
  three regions; 129 BioProjects; China ST34 mcr-1 counts
  (351 Human, 125 Human Hybrid, 13 Chicken_related, 4 Chicken_related Hybrid).

Usage

  Run from a working directory containing both files:

      Rscript R_regional_sensitivity.R

  The script creates R_sensitivity_output/.

Outputs

  regional_counts.csv
      Regional sample sizes, Hybrid counts, pooled weights, and
      source-specific Hybrid rates for Asia, North America, and Europe.

  leave_one_BioProject_out.csv
      Leave-one-BioProject-out estimates for each of the 129 BioProjects.

  analysis_summary.csv
      Summary of all-isolate rates, standardized rates and differences,
      bootstrap intervals, leave-one-project-out range, valid bootstrap
      replicates, and China ST34 mcr-1 counts.

  R_session_info.txt
      R session information.

Analysis design

  Overall cohort
      1,514 genomes: 879 Human, 635 Chicken_related, 223 Hybrid.

  Three-region subset
      Asia, North America, and Europe; 1,442 genomes; 202 Hybrid.

  Geographic standardization
      Fixed pooled regional weights from the three-region subset applied
      to source-specific Hybrid rates.

  Alternative weights
      Equal three-region weights; Asia and North America only.

  BioProject cluster bootstrap
      30,000 resamples of 129 BioProjects with replacement; resampled
      copies retain all source-region counts of the selected BioProject;
      empty source-region replicates are omitted; seed 20261013.

  Independent-isolate parametric binomial bootstrap
      Fixed source-region sample sizes; seed 20261014.

  Leave-one-BioProject-out
      Each BioProject removed once; original pooled regional weights fixed.

  China ST34 mcr-1 subset
      Counts for Country = China, MLST = 34, mcr subtype = mcr-1.