Complete mcr-positive plasmid analysis: code and data

Contents

  00_admin/filter_TADB2_component_hits.py
      TADB2 component filtering. Reads the raw BLASTP table, applies the
      strict filter (pident > 80, qcovhsp > 80, E-value < 1e-5), and joins
      target metadata.

  00_admin/validate_TADB2_typeII_modules_fixed.py
      TADB2 physical-pair validation. Reads the full passing-hit table,
      requires distinct proteins with compatible toxin/antitoxin assignments
      sharing a TADB reference pair ID, and applies the adjacency criterion
      (no intervening translated CDS and maximum 500 bp genomic gap,
      including circular-origin adjacency). Candidate modules that reuse a
      protein are resolved greedily by combined evidence score.

  scripts/01_plasmid_summary.py
      Per-plasmid species, ARG, and backbone summaries from existing
      937-plasmid annotation tables.

  scripts/02_TA_ARG_associations.py
      TA-ARG Spearman correlations and backbone strata from the same
      937-plasmid annotation tables.

  scripts/03_TADB2_blastp.sh
      makeblastdb and BLASTP commands for TADB2 query proteins against
      predicted proteins from 937 complete plasmids.

  scripts/04_mcr_candidate_blastn.sh
      Relaxed BLASTN candidate search for mcr sequences.

  scripts/05_resfinder_confirmed_smoketest.sh
      ResFinder database build and three-plasmid ABRicate smoke test.

  data/plasmid_metadata_937.tsv
      937-plasmid metadata.

  data/plasmid_non_mcr_ARG_937.tsv
      Non-mcr ARG counts and combinations.

  data/plasmid_TA_937.tsv
      TA module annotations.

  SOURCE_SCRIPT_MANIFEST.tsv
      Inventory of scripts and prerequisite data.

  RUN_ORDER_TADB2.txt
      Execution order for the TADB2 scripts.

  requirements.txt
      Python dependencies.

Requirements

  Python 3.8 or later
  numpy, pandas, scipy (see requirements.txt)
  For shell scripts: makeblastdb, blastp, blastn, abricate

Data files

  data/plasmid_metadata_937.tsv contains 937 plasmid accessions.
  data/plasmid_non_mcr_ARG_937.tsv contains non-mcr ARG annotations for the
  same accessions.
  data/plasmid_TA_937.tsv contains TA module annotations for the same
  accessions.

  The 937-plasmid RefSeq collection is a multi-species collection. It is not
  one-to-one matched to the 1,514 Salmonella isolates analyzed elsewhere.

Run order for TADB2 scripts

  The TADB2 scripts are run from a project root containing the 06_TA_systems/
  input files. The original working directory structure used 00_admin/ and
  06_TA_systems/. The original private working directory has been removed;
  run scripts in a copy of the project.

  1. bash scripts/03_TADB2_blastp.sh
  2. python3 00_admin/filter_TADB2_component_hits.py
  3. python3 00_admin/validate_TADB2_typeII_modules_fixed.py

  The reference preparation step originally used
  python3 00_admin/prepare_TADB2_typeII_reference.py. This source file is not
  included; prepared reference FASTA and metadata are prerequisites.
  The final TA family-curation step originally used
  python3 00_admin/curate_TADB2_selected_pair_families.py. This source file
  is not included.

Downstream statistics

  python3 -m pip install -r requirements.txt
  python3 scripts/01_plasmid_summary.py --data data --out results
  python3 scripts/02_TA_ARG_associations.py --data data --out results

Software and parameters

  TADB 2.0 Type II BLASTP
    Query: prepared toxin and antitoxin reference proteins from TADB 2.0
    Subject database: proteins predicted from 937 complete plasmids
    BLASTP: -evalue 1e-5 -seg no -max_hsps 1 -max_target_seqs 5000
            -num_threads 48
    Output fields include qcovhsp
    Filtering: pident > 80, qcovhsp > 80, E-value < 1e-5
      (implemented in 00_admin/filter_TADB2_component_hits.py)
    Pair validation: distinct proteins, same plasmid, no intervening
      translated CDS, inter-CDS spacing <= 500 bp (overlaps allowed),
      >=1 shared reference TADB pair ID
      (implemented in 00_admin/validate_TADB2_typeII_modules_fixed.py)
    Reference deduplication: within a component and pair ID, identical amino
      acid sequences are collapsed; different protein variants are retained.
    Historical intermediate counts: 538 plasmids with >=1 homologous Type II
      protein; 583 preliminary candidate pairs; 31 distinct reference pair IDs
      in the stricter set. These counts use different denominators.

  mcr BLASTN screen
    Query: de-duplicated primary mcr sequence set compiled from a ResFinder
      colistin reference snapshot (109 entries before sequence processing).
    BLASTN candidate retrieval: -perc_identity 85 -evalue 1e-20
    Final inclusion threshold: >=90% nucleotide identity, >=80% query
      coverage, accession deduplication, and subtype checks.
    The final filter is a separate downstream step. Its original Python
      source is not included in this package.

  ResFinder and replicon annotation
    ResFinder: ABRicate against frozen local database resfinder_eecf0aa20759,
      minimum identity and coverage both 80%.
    scripts/05_resfinder_confirmed_smoketest.sh is a three-plasmid smoke
      test, not a full 937-plasmid batch command.
    MOB-typer v3.1.9 was used for replicon, relaxase, oriT, and MPF typing.
    A separate PlasmidFinder executable or script invocation is not included
      in this package.

Outputs

  scripts/01_plasmid_summary.py writes:
    plasmid_host_species_summary.tsv
    selected_co_resistance_plasmids.tsv
    plasmid_backbone_summary.tsv
    summary_counts.txt

  scripts/02_TA_ARG_associations.py writes:
    TA_ARG_overall_spearman.tsv
    TA_ARG_by_backbone.tsv

  TADB2 scripts write outputs under 06_TA_systems/03_blast/results/ and
  06_TA_systems/04_pair_validation/.

Testing status

  The Python scripts were syntax-checked. They were not rerun on the original
  raw protein and BLAST files in this package.

Manifest

  SOURCE_SCRIPT_MANIFEST.tsv lists each analysis step, file or artifact,
  status, and purpose.