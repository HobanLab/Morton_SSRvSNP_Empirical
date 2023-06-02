#!/bin/bash

# %%% OUTLINE %%%

# This script runs the populations module for the QUAC optimized de novo assembly (m7, M4, n4, gt-alpha 0.01).
# The populations module filters loci according to missing data/allele frequencies, and generate different file outputs
# which are needed for downstream analyses. It also allows users to specify which samples to include in final outputs
# by listing those samples (and their respective populations) in the popmap file.

# Multiple population module scripts are called: the outputs of these are used for different aspects of downstream analyses
# (for instance, ex situ conservation, population clustering, etc.). The parameters that are varied are listed below:
#
# + Missing data filters: R0 (no data filter, and no flag used) and R80 (loci shared amongst 80% of samples)
# + Wild: these scripts analyze only the wild samples, and utilize the QUAC_popmap_wild population map
# + gardenProv: these scrips are specific to QUAC, and are used for assessing sample garden provenances
# + Number of SNPs per locus: we specify writing only the first SNP per locus (1SNP), or writing SNPs haplotype-wise (HapSNPs)

# In general, we calculate divergence from HWE for each locus, and specify creating most possible
# filetype outputs for each population module run. We fix the number of threads at 16.

# Specify complete path to the results of the denovo_map.pl run. Population module runs are based off of the outputs from this run
QUAC_DN_DIR="/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/originalOutput"

# %%% POPULATIONS MODULE CALLS %%%

# GARDEN AND WILD %%%
# Complete
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap -t 16 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap -t 16 -R 80 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all

# HAPLOTYPE-WISE SNPS %%%
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap -t 16 -hwe -H --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap -t 16 -R 80 -hwe -H --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all

# WILD ONLY %%%
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap_Wild -t 16 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap_Wild -t 16 -R 80 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all

# GARDEN PROVENANCE %%%
# The outputs from these runs are used to assess proper assignment of garden samples back to wild source populations (using STRUCTURE)
# The sample sets used for these samples are Subset garden and wild samples, with garden samples grouped according to their documented source population
populations -P $QUAC_DN_DIR -O ./ -M ./QUAC_popmap_gardenProv -t 16 -R 80 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
