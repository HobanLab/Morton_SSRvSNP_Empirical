## RADseq Analyses

This folder contains the scripts and files used to generate _de novo_ assemblies and reference alignments of 
the NextRAD data, including the scripts used to run relevant Stacks commands. Four folders are included.

### filtering_QC
This folder contains scripts used to run the [process_radtags](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php)
command, as well as analysis of reads (before and after cleaning with process_radtags) using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### denovoAssembly
This folder contains the BASH and R scripts used to run parameter optimization for the [Stacks denvo_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php)
script. Then, using optimized parameters (for both QUAC and QUBO), _de novo_ assemblies were built. 

### referenceAlignment
This folder contains the BASH scripts used to download and index the reference genomes used in our analyses, as well as
the scripts used to align our NextRAD data to those reference genomes. The scripts calling the [Stacks ref_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php)
commands are also included in this folder, as well as scripts generating outputs that report the quality of reference alignment.

### hybridAnalysis
This folder contains the BASH scripts used for running the hybrid analysis, in which contigs assembled in optimized
_de novo_ assemblies were aligned to reference genomes. We ultimately did not utilize these outputs (we only compared _de novo_
assembly reads and reference alignments), but the scripts are kept here for documentation purposes.

### duplicateAnalysis.R
This R script measures the number of loci lost or gained when duplicate samples are removed from either _de novo_ or
reference aligned datasets. The quantification of duplicate loci was used to assess the quality of our different
SNP processing approaches (since the number of loci identified between duplicate samples is indicative of sequencing error). 
