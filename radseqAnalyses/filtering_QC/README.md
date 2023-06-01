## Filtering and QC Steps 

This folder contains the BASH scripts used to run the FastQC software, for analyzing reads. It also includes the 
process_radtags commands used to clean NextRAD data, prior to downstream analyses (de novo assembly or reference
alignment). 

### fastQC
This folder contains the output of the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) software. 
The command used to run this software is provided in the `FastQC_run.sh` script. 

The folder contains 3 subdirectories: one for the FastQC outputs of our raw NextRAD data ("Prefiltering"),
and two from our data cleaned with the Stacks `process_radtags` command. We ran `process_radtags`
*without* adapter trimming ("Untrimmed") and also *with* adapter trimming ("Trimmed"). More information 
is included in the [processRADtags.sh](https://github.com/akoontz11/Morton_SSRvSNP_Empirical/blob/main/radseqAnalyses/filtering_QC/dataCleaning.sh) script.

### processRADtags.sh
This script is a collection of various commands used to clean our NextRAD data. These range from pulling our read
files from the sequencing center sharesite, to changing file limits on our server to run `process_radtags`, to the
`process_radtags` commands themselves. 

The last of these is arguably the most important. For our study, we ultimately used the `process_radtags` command 
_with_ trimming, for all downstream analyses (de novo assemblies/reference alignments).

### analyzeReads.R
This script analyzes the outputs of the `stacks-dist-extract process_radtags.log per_barcode_raw_read_counts` command
(that is, it assesses the quality of the outputs from process_radtags).
Specifically, it calculates means/standard deviations of raw read counts.

### Barcodes.txt
The Barcodes.txt file was used by the `process_radtags` commands to sort samples by their species
(QUBO = _Quercus boyntonii_; QUAC = _Quercus acerifolia_)and collection origin (Garden, G, or Wild, W).
Samples were originally named according the the combination of adatpers used to identify them; 
this file renames them according to their IMLS/SH_Q database names (databases internal to the Hoban Lab).
