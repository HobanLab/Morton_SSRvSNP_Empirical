#!/bin/bash

# %%% OUTLINE %%%

# This script runs the populations module for the QUBO reference alignment (with the Q. robur reference genome).
# The populations module filters loci according to missing data/allele frequencies, and generate different file outputs
# which are needed for downstream analyses. It also allows users to specify which samples to include in final outputs
# by listing those samples (and their respective populations) in the popmap file.

# Multiple population module scripts are called: the outputs of these are used for different aspects of downstream analyses
# (for instance, ex situ conservation, population clustering, etc.). The parameters that are varied are listed below:
#
# + Missing data filters: R0 (no data filter, and no flag used) and R80 (loci shared amongst 80% of samples)
# + Wild: these scripts analyze only the wild samples, and utilize the QUBO_popmap_wild population map
# + Subset: these scripts analyze only the samples shared between MSAT and SNP datasets (Subset samples).
#           Because there are differences in the wild samples used between MSAT and SNP datasets (unlike QUAC),
#           we need to utilize a popmap file which only contains the wild samples present in the SNP datasets.
# + Number of SNPs per locus: we specify writing only the first SNP per locus (1SNP), or writing SNPs haplotype-wise (HapSNPs)

# In general, we calculate divergence from HWE for each locus, and specify creating most possible
# filetype outputs for each population module run. We fix the number of threads at 16.

# Specify complete path to the results of the denovo_map.pl run. Population module runs are based off of the outputs from this run
QUBO_REF_DIR="/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/originalOutput"

# %%% POPULATIONS MODULE CALLS %%%

# GARDEN AND WILD %%%
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap -t 16 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap -t 16 -R 80 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all

# HAPLOTYPE-WISE SNPS %%%
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap -t 16 -hwe -H --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap -t 16 -R 80 -hwe -H --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all

# WILD ONLY %%%
# Complete
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap_Wild -t 16 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap_Wild -t 16 -R 80 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
# Subset
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap_Wild_Subset -t 16 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all
populations -P $QUBO_REF_DIR -O ./ -M ./QUBO_popmap_Wild_Subset -t 16 -R 80 -hwe --write-single-snp --fstats --fasta-samples --vcf --genepop --structure --plink --phylip-var-all

