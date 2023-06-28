# Overview
This repository contains the code used in the SSR vs. SNP marker comparison project,
one of the initiatives in the [GCCO](https://www.globalconservationconsortia.org/gcc/oak/) group. 
We compare microsatellites and SNPs from a RADseq analysis, and calculate how each marker leads to 
different measurements of _ex situ_ representation and population clustering.
We do this using two different oak species, _Quercus boyntonii_ (or QUBO; Boynton oak, in Alabama, in the white oak
section of the genus) and _Quercus acerifolia_ (or QUAC; a red oak, section Lobatae). 
Both of these species are included on the IUCN Red List.

# Repository layout
Below we describe the general structure of the repository and its directories; more information can be found in each
folder's respective READMEs. Generally, this repo will contain the scripts used to generate data for different steps in our analyses
and any other input files required for those scripts to run, such as parameter values, sample lists, 
and input files. However, many input files cannot be included in this repository due to GitHub file size limits
(for instance, no genind files, and no input or output STRUCTURE files).

The [radAnalysis](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/radAnalysis) folder represents the first part of this process, which involves QC and SNP calling (via _de novo_
assembly and reference alignment). The [exSituRep](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/exSituRep) folder contains the scripts for calculating _ex situ_ representation and
resampling analyses, and the [popAnalysis](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/popAnalysis) folder contains the scripts for clustering and calculating population level data.

Often, subfolders will be named QUAC and QUBO: these refer to the files used for the different 
species being analyzed. Files will often match one another for both species (but not always!)

## radAnalysis
The files in this folder describe the steps used to generate SNP datasets for downstream _ex situ_
representation analyses. Assemblies and loci derived from reference alignments were built using the
[Stacks](https://catchenlab.life.illinois.edu/stacks/) software.

### denovo
This contains the scripts and parameter files used to explore the parameter space of a Stacks 
_de novo_ assembly. It contains 2 folders: one for parameter optimization ([paramOpt](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/radAnalyses/denovo/paramOpt)), 
and one for SNP calling ([finalAssemblies](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/radAnalyses/denovo/finalAssemblies)).

### filtering_QC
This folder contains various scripts used to clean the raw NextRAD data (generated from the
sequencing company SNPsaurus). It also includes the FastQC scripts used to assess data quality.

### reference
This folder contains the scripts used to download and index reference genomes, build reference alignments, and call SNPs from those alignments. 
It contains 2 folders: one for the scripts used in building reference alignments ([refAlign](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/radAnalyses/reference/refAlignment)), 
and one for SNP calling ([processLoci](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/tree/main/radAnalyses/reference/processLoci)).

## exSituRep
This folder contains the scripts used to run two analyses concerned with the _ex situ_ representation of QUAC/QUBO, described below.

### exSituRepresentation
The scripts in this folder measure how well garden samples represent the allelic diversity of wild populations, in different contexts. 

There are 2 scripts for _ex situ_ representation calculations. The first, exSituRepresentation.R, is the primary script, the results of which are presented in our
final reports. The second, exSituRepresentation_Expanded.R, is more exploratory: it examines the impacts of different filters (mostly from the Stacks
[populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) module) on _ex situ_ representation.

The functions used for both scripts, as well as for resampling analyses, are declared in the [functions_exSituRepresentation.R](https://github.com/HobanLab/Morton_SSRvSNP_Empirical/blob/main/exSituRep/functions_exSituRepresentation.R) script.
These functions are nested, but ultimately allow for the calculation of _ex situ_ representation given a single genind file with two populations: "wild" and "garden".
The functions used for resampling are also nested, in order to allow for parallel processing of large SNP datasets.

### Resampling
The scripts in this folder run the resampling analyses, in which wild samples are iteratively selected randomly, allelic representation
of those randomly selected samples is calculated, and this process is repeated to determine a minimum sampling size estimate for different 
thresholds of required wild allelic diversity. This folder also contains a script for generating images used in the manuscript for this project.

## popAnalysis
This folder contains R scripts used for analyzing QUAC/QUBO populations using different approaches common in 
population genetic analyses. Scripts in here are used to calculate statistical measures (Fst, heterozygosity, allelic richness)
and run clustering analyses (STRUCTURE and DAPC). Note that the Fst calculations are separated from the He and Ar calculations.

### STRUCTURE
This folder contains the script used to create STRUCTURE plots (one for the project in general, ane one for the manuscript).
It also contains a script that's used to convert the microsatellite genind files to STRUCTURE. This script uses the function [genind2structure](https://zenodo.org/record/846816),
which was originally written by [Lindsay Clark](https://github.com/lvclark).

# Data and Contact
A link to the raw data for this analysis can be found on the [NCBI SRA](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)

For questions about this dataset or the scripts included here, contact Austin Koontz.

