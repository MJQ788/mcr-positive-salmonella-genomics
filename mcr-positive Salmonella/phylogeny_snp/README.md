Core-SNP phylogeny and <=10-SNP network analysis



Overview



&#x20; This package contains scripts for core-genome SNP phylogeny and

&#x20; <=10-SNP network construction for 724 Typhimurium-complex genomes.

&#x20; Inputs are assembled genome contigs. Outputs include a core

&#x20; whole-reference alignment, a recombination-filtered SNP alignment, a

&#x20; maximum-likelihood tree, a pairwise SNP distance matrix, and connected

&#x20; components defined by <=10-SNP edges.



Contents



&#x20; 01\_snippy\_and\_core.sh

&#x20;     Runs Snippy on assembled contigs and combines isolate outputs with

&#x20;     Snippy-core.



&#x20; 02\_gubbins\_snp\_sites.sh

&#x20;     Runs snippy-clean\_full\_aln, Gubbins, and snp-sites -c.



&#x20; 03\_iqtree.sh

&#x20;     Runs IQ-TREE maximum-likelihood inference with K3Pu+ASC+R2 and

&#x20;     1,000 ultrafast bootstrap replicates.



&#x20; 04\_snp\_dists.sh

&#x20;     Calculates pairwise SNP distances with SNP-dists.



&#x20; 05\_network\_connected\_components.py

&#x20;     Builds <=10-SNP undirected edges and connected components from a

&#x20;     SNP-dists matrix.



Software and parameters



&#x20; Snippy with --ctgs; 50 CPUs.

&#x20; Snippy-core with a single reference genome.

&#x20; Gubbins with 40 threads.

&#x20; SNP-sites with -c, retaining SNP columns containing only A/C/G/T.

&#x20; IQ-TREE model K3Pu+ASC+R2; 1,000 ultrafast bootstrap replicates;

&#x20;   3 threads.

&#x20; SNP-dists on clean.core.aln.

&#x20; <=10-SNP network: an edge joins two isolates when their core-SNP

&#x20;   distance is <=10; connected components are defined by graph

&#x20;   connectivity, not by requiring all pairwise distances within a

&#x20;   component to be <=10.



Input layout



&#x20; Run from the analysis project directory. Example:



&#x20; export INPUT\_DIR=genomes

&#x20; export OUTPUT\_DIR=snippy\_output

&#x20; export REFERENCE=reference/re.gbk

&#x20; export SNIPPY\_CPUS=50

&#x20; export GUBBINS\_THREADS=40



&#x20; INPUT\_DIR contains assembled .fna files.

&#x20; REFERENCE is the reference genome file.

&#x20; Optional metadata for network annotation is a CSV or TSV with isolate

&#x20; ID and source columns.



Usage



&#x20; bash scripts/01\_snippy\_and\_core.sh

&#x20; bash scripts/02\_gubbins\_snp\_sites.sh core.full.aln clean.full.aln gubbins\_result clean.core.aln

&#x20; bash scripts/03\_iqtree.sh clean.core.aln iqtree\_final\_tree

&#x20; bash scripts/04\_snp\_dists.sh clean.core.aln distances.tsv

&#x20; python3 scripts/05\_network\_connected\_components.py --distances distances.tsv --metadata metadata.tsv --id-column Assembly --source-column "Isolation source"



Outputs



&#x20; core.full.aln

&#x20;     Core whole-reference alignment from Snippy-core.



&#x20; clean.full.aln

&#x20;     Snippy-cleaned full alignment.



&#x20; gubbins\_result.filtered\_polymorphic\_sites.fasta

&#x20;     Recombination-filtered polymorphic sites from Gubbins.



&#x20; clean.core.aln

&#x20;     SNP alignment after snp-sites -c.



&#x20; iqtree\_final\_tree.\*

&#x20;     IQ-TREE output files.



&#x20; distances.tsv

&#x20;     Pairwise SNP distance matrix.



&#x20; snp\_network/edges\_leq10.tsv

&#x20;     Undirected edges with distance <=10.



&#x20; snp\_network/components.tsv

&#x20;     Connected components and source counts.



&#x20; snp\_network/mixed\_source\_nodes.tsv

&#x20;     Nodes in components containing both human and chicken-related

&#x20;     sources.



Software documentation



&#x20; Snippy: https://github.com/tseemann/snippy

&#x20; SNP-sites: https://github.com/sanger-pathogens/snp-sites

&#x20; Gubbins: https://github.com/nickjcroucher/gubbins

&#x20; IQ-TREE: https://www.iqtree.org/doc/Substitution-Models

&#x20; SNP-dists: https://github.com/tseemann/snp-dists

