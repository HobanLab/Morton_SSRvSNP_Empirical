## Filtering and QC Steps 

This folder contains the BASH scripts and commands used to clean the NextRAD data prior to any Stacks commands or alignment to 
reference genomes. 

### fastQC
This folder contains the output of the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) software. 
The command used to run this software is provided in the `FastQC_run.sh` script. 

The folder contains 3 subdirectories: one for the FastQC outputs of our raw NextRAD data, 
and two from our data cleaned with the Stacks `process_radtags` command. We ran `process_radtags`
*without* adapter trimming ("Untrimmed") and also *with* adapter trimming ("Trimmed"). More information 
is included in the [dataCleaning.sh](https://github.com/akoontz11/Morton_SSRvSNP_Empirical/blob/main/radseqAnalyses/filtering_QC/dataCleaning.sh) script.

### preliminarySubset
This folder contains the scripts and Stacks files used to analyze a 12 sample subset
of QUBO and QUAC samples together. The purpose of the preliminarySubset analysis was to check for
basic assumptions of our data: making sure replicates matched, and making sure samples within the same
populations were more similar to each other than samples across populations. 

There are subfolders for the scripts passed to `stacks`, STRUCTURE (for analyzing the Stacks output). The `analysis`
folder  contains R scripts used for analyzing the Stacks outputs.

### dataCleaning.sh
This script is a collection of various commands used to clean our NextRAD data. These range from pulling our read
files from the sequencing center sharesite, to changing file limits on our server to run `process_radtags`, to the
`process_radtags` commands themselves. 

The last of these is arguably the most important. For our study, we ultimately used the `process_radtags` command 
_with_ trimming, for all downstream analyses (de novo assemblies/reference alignments).

### analyzeReads.R
This script analyzes the outputs of the `stacks-dist-extract process_radtags.log per_barcode_raw_read_counts` command.
Specifically, it calculates means/standard deviations of raw read counts.

### Barcodes.txt
The Barcodes.txt file was used by the `process_radtags` commands to sort samples by their species
(QUBO = _Quercus boyntonii_; QUAC = _Quercus acerifolia_)and collection origin (Garden, G, or Wild, W).
Samples were originally named according the the combination of adatpers used to identify them; 
this file renames them according to their IMLS/SH_Q database names (databases internal to the Hoban Lab).
