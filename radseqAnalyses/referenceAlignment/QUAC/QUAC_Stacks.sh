#!/bin/bash

# This script contains the Stacks ref_map.pl and populations commands for filtered, reference aligned QUAC reads

# %%% IDENTIFYING REFERENCE ALIGNMENT LOCI %%%
# The only modules that run for ref_map.pl are gstacks and populations;
# we pass the --rm-pcr-duplicates and --gt-alpha arguments to gstacks in order to mirror the parameters used for de novo assembly.
# Paired end reads were passed to the GMAP software, so Stacks should automatically detect that it's working with paired end data

# Command for processing reference alignments, specifying the input (--samples) and output directories
ref_map.pl -T 28 --gt-alpha 0.01 --rm-pcr-duplicates --samples /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/REF_QURU/GSNAP4 --popmap ./QUAC_popmap -o ./output/

# After ref_map.pl script has completed, analyze the gstacks outputs to assess reference alignment quality
stacks-dist-extract gstacks.log.distribs bam_stats_per_sample >> bamStats.csv

# %%% PROCESSING STACKS OUTPUTS %%%
# Each line below includes a different command sent to the Stacks populations module, for processing the output of the ref_map.pl command above
# The first line specifies an R paramter of 0 (no filter on missing data). The second line specifies an R parameter of 80 (only loci seen within
# 80% of samples or more will be included in final dataset). This commands are run from folders within the output directory of the above
# ref_map.pl command, which is why the -P parameter specifies the directory above the run directory.

# Other arguments are as follows:
# 15 threads; no minimum allele frequency is specified; only the 1st SNP per locus is written; calculate divergence from HWE for each locus
# Output files: Haplotypes in FASTA format; VCF files; GenePop files; STRUCTURE files; PLINK files; PHYLIP files

# No missing data filter
populations -P ../ -O ./ -M ./QUAC_popmap_GardenWild -t 15 --hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
# R80 missing data fileter
populations -P ../ -O ./ -M ./QUAC_popmap_GardenWild -t 15 -R 80 --hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
