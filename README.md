# Overview
This repository contains the code used in the SSR vs. SNP marker comparison project,
one of the initiatives in the [GCCO](https://www.globalconservationconsortia.org/gcc/oak/) group. 
We compare microsatellite and SNPs from a RADseq analysis, and calculate how each marker leads to 
different measurements of ex situ representation.
We do this using two different oak species, _Quercus boyntonii_ (or QUBO; Boynton oak, in Alabama, in the white oak
section of the genus) and _Quercus acerifolia_ (or QUAC; a red oak, section Lobatae). 
Both of these species are included on the IUCN Red List.

Below we describe the general structure of the repository and it's directories. Generally,
this repo will contain the scripts used to generate data for different steps in our analyses
and any other input files required for those scripts to run, such as parameter values, sample lists, 
and input files. However, many input files cannot be included in this repository due to GitHub file size limits
(for instance, there are no genind files, and no input or output STRUCTURE files).

Often, subfolders will be named QUAC and QUBO: these refer to the files used for the different 
species being analyzed. Files will often match one another for both species (but not always!)

## radseqAnalyses
The files in this folder describe the steps used to generate SNP datasets for downstream _ex situ_
representation analyses. Assemblies and loci derived from reference alignments were built using the
[Stacks](https://catchenlab.life.illinois.edu/stacks/) software. More information can be found in the folder's 
specific README, but a rough outline is provided below.

### denovo
This contains the scripts and parameter files used to explore the parameter space of a Stacks 
_de novo_ assembly. It contains 2 folders: one for parameter optimization ([paramOpt]()), and one for SNP calling ([finalAssemblies]()).

More information can be foud in the README for this folder.

### filtering_QC
This folder contains various scripts used to clean the raw NextRAD data (generated from the
sequencing company SNPsaurus). The fastQC folder contains the output of the `fastqc` program,
which was run on our raw data ("Prefiltering"), our data cleaned with `process_radtags` *without* adapter trimming
("Untrimmed"), and our data cleaned with `process_radtags` *with* adapter trimming. 

The Barcodes.txt file was used by the `process_radtags` command to sort samples by their species
and collection origin (Garden, G, or Wild, W). The dataCleaning script contains the command `process_radtags`
commands sent to Stacks, as well as other commands used to generate and organize Stacks output.

More information can be foud in the README for this folder.

### reference
This folder contains the scripts used to download and index reference genomes, build reference alignments, and call SNPs from those alignments. 
It contains two folders: one for the scripts used in building reference alignments ([refAlign]()), and one for SNP calling ([processLoci]()).

More information can be foud in the README for this folder.

## popAnalyses
This folder contains R scripts used for analyzing QUAC/QUBO populations using different approaches common in 
population genetic analyses. Scripts in here are used to calculate statistical measures (Fst, heterozygosity, allelic richness)
and run clustering analyses (STRUCTURE and DAPC).

## exSituRepresentation
This folder contains the scripts used to run two analyses concerned with the _ex situ_ representation of QUAC/QUBO. 
These are described below. The functions_exSituRepresentation.R script declares functions used in both of these analyses.

### exSituRepresentation
The scripts in this folder measure how well garden samples represent the allelic diversity of wild populations, in different contexts. 

The exSituRepresentation_SNPs.R script measures wild allelic representation in gardens for the QUAC/QUBO SNP datasets, but explores **many**
different filters on those datasets. 

The exSituRepresentation_MSATs.R script measures wild allelic representation in gardens for the QUAC/QUBO SNP *AND* MSAT datasets,
both complete (all samples for each dataset) and subset (only samples shared across SNP and MSAT studies, for each species). 

The exSituRepresentation_QUHA.R script measures wild allelic representation in the same datasets as the exSituRepresentation_MSATs.R script.
It's meant to serve as a validation for the values generated in that script.

### Resampling
The scripts in this folder run the resampling analyses, in which wild samples are iteratively selected randomly, allelic representation
of those randomly selected samples is calculated, and this process is repeated to determine a minimum sampling size for different 
thresholds of required wild allelic diversity.

## Contact
For questions about this dataset, open an Issue or contact Austin. 

A link to the raw data for this analysis can be found on the [NCBI SRA](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)

