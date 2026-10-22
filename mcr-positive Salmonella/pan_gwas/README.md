Pan-GWAS analysis scripts for mcr-positive Salmonella



Contents



Core scripts



&#x20; 00\_generate\_filtered\_rtabs.py

&#x20;     Generates conservative and strict feature-filtered Rtab matrices

&#x20;     from the full Rtab and Panaroo outputs.



&#x20; 01\_pyseer\_workflow.sh

&#x20;     Runs Pyseer for the full, conservative, and strict feature sets.



&#x20; 02\_full\_gwas\_qc\_fdr.py

&#x20;     Applies main-QC filtering and BH-FDR correction to the full Pyseer

&#x20;     output.



&#x20; 03\_stage1.py

&#x20;     Calculates feature frequencies and selects exact co-presence block

&#x20;     representatives.



&#x20; 05\_candidate\_prioritization.py

&#x20;     Applies annotation-aware prioritization (v2) to Stage-1 top-block

&#x20;     representatives.



Additional scripts



&#x20; 04\_prepare\_annotated\_topblocks.py

&#x20;     Applies the >=20 percentage-point / >=30% Hybrid prevalence screen

&#x20;     and merges Panaroo annotations into Stage-1 top-block

&#x20;     representatives.



&#x20; 06\_filtered\_gwas\_qc\_fdr.py

&#x20;     Applies main-QC filtering and per-analysis BH-FDR correction to the

&#x20;     conservative and strict Pyseer outputs.



Input layout



&#x20; Run all commands from the project root. Input files are arranged as:



&#x20; panGWAS\_input/

&#x20;   pyseer\_Hybrid\_phenotypes.tsv

&#x20;   pyseer\_coreSNP\_distances.tsv

&#x20;   pyseer\_gene\_presence\_absence\_min1pct.Rtab

&#x20;   panGWAS\_full\_pyseer\_alltests.tsv

&#x20; panaroo\_result\_moderate/

&#x20;   gene\_presence\_absence.csv

&#x20;   gene\_data.csv

&#x20; pyseer\_results/

&#x20; scripts/



Software and parameters



&#x20; Pyseer v1.4.2; fixed-effects association of gene presence-absence with

&#x20; Hybrid status; population correction using four MDS coordinates from

&#x20; core-SNP distances; --min-af 0.01 --max-af 0.99; --print-samples.

&#x20; --cpu controls parallelization only.



&#x20; Panaroo v1.5.2 (moderate cleaning) and Prokka v1.14.6 upstream.



&#x20; Python 3 with pandas and numpy.



&#x20; A row enters BH-FDR only if its LRT P value is finite and its notes

&#x20; field is empty or none.



&#x20; BH correction is applied independently to each feature set.



&#x20; Conservative filtering excludes annotation-flagged AMR/efflux and VF

&#x20; features, plus exact co-presence blocks containing those seed features.



&#x20; Strict filtering additionally excludes replication- and

&#x20; methylation-associated features and their exact co-presence blocks.



&#x20; Filtering keywords include \\bant and \\bsip.



Command sequence



&#x20; # 1. Generate the two feature-filtered matrices from Panaroo outputs.

&#x20; python3 scripts/00\_generate\_filtered\_rtabs.py --project-root .



&#x20; # 2. Run the three Pyseer analyses.

&#x20; bash scripts/01\_pyseer\_workflow.sh full

&#x20; bash scripts/01\_pyseer\_workflow.sh conservative

&#x20; bash scripts/01\_pyseer\_workflow.sh strict



&#x20; # 3. Main-QC and BH-FDR of the full Pyseer output.

&#x20; (cd panGWAS\_input \&\& python3 ../scripts/02\_full\_gwas\_qc\_fdr.py)



&#x20; # 4. Stage-1 frequency calculations and co-presence block

&#x20; #    representatives.

&#x20; python3 scripts/03\_stage1.py --project-root .



&#x20; # 5. Top-block screening and Panaroo annotation merge.

&#x20; python3 scripts/04\_prepare\_annotated\_topblocks.py --project-root .



&#x20; # 6. Annotation-aware prioritization (v2).

&#x20; python3 scripts/05\_candidate\_prioritization.py



&#x20; # 7. Main-QC and BH-FDR of the conservative and strict Pyseer

&#x20; #    outputs.

&#x20; python3 scripts/06\_filtered\_gwas\_qc\_fdr.py --project-root .



Outputs



&#x20; panGWAS\_input/pyseer\_gene\_presence\_absence\_min1pct.filtered\_conservative.Rtab

&#x20;     Conservative feature-filtered Rtab matrix.



&#x20; panGWAS\_input/pyseer\_gene\_presence\_absence\_min1pct.filtered\_strict.Rtab

&#x20;     Strict feature-filtered Rtab matrix.



&#x20; panGWAS\_input/panGWAS\_full\_pyseer\_mainQC\_notesNone\_with\_FDR.tsv

&#x20;     Main-QC-passed full Pyseer table with BH-FDR q values.



&#x20; panGWAS\_input/panGWAS\_full\_pyseer\_mainQC\_FDR0.05.tsv

&#x20;     Full Pyseer features with q < 0.05.



&#x20; pyseer\_results/Candidate\_Prioritization\_Stage1\_All.tsv

&#x20;     Stage-1 features with frequencies and block assignments.



&#x20; pyseer\_results/Candidate\_Prioritization\_Stage1\_Hybrid\_Enriched.tsv

&#x20;     Hybrid-enriched Stage-1 features.



&#x20; pyseer\_results/Candidate\_CoPresence\_Block\_Summary.tsv

&#x20;     Co-presence block summary.



&#x20; pyseer\_results/Candidate\_Prioritization\_Stage1\_Hybrid\_Enriched\_TopBlocks\_Annotated.tsv

&#x20;     Top-block representatives with Panaroo annotations.



&#x20; pyseer\_results/Candidate\_Prioritization\_Stage1\_Hybrid\_Enriched\_TopBlocks\_Annotated\_Reclassified\_v2.tsv

&#x20;     Annotation-aware prioritization output.



&#x20; pyseer\_results/Candidate15\_Initial\_Selection\_Audit.tsv

&#x20;     Initial 15 candidate selection audit table.



&#x20; pyseer\_results/Candidate15\_Initial\_Gene\_IDs.txt

&#x20;     Initial 15 candidate gene IDs.



&#x20; pyseer\_results/filtered\_nondefinition\_gwas/panGWAS\_filtered\_conservative\_mainQC\_FDR.tsv

&#x20;     Conservative-filtered Pyseer output with main-QC and BH-FDR.



&#x20; pyseer\_results/filtered\_nondefinition\_gwas/panGWAS\_filtered\_strict\_mainQC\_FDR.tsv

&#x20;     Strict-filtered Pyseer output with main-QC and BH-FDR.



Software documentation



&#x20; Pyseer: https://github.com/mgalardini/pyseer

&#x20; Panaroo: https://github.com/gtonkinhill/panaroo

&#x20; Prokka: https://github.com/tseemann/prokka

