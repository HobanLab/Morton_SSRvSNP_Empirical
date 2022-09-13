#!/bin/bash

# %%%%% PRELIMINARY ANALYSIS %%%%%
# We assembled loci de novo using a representative subset of all of our samples (both species, both garden and wild samples). This was to ensure some basic checks about our data:
# 1. Replicates should match one another
# 2. Samples from the same populations should resemble one another more than they do samples from other populations (for wild samples).

# Run standalone Stacks populations program, and generate VCF, STRUCTURE, genepop, and Fst outputs. Use all loci, but only analyze first SNP per locus
populations -P ./output/ -O ./output/populations_all/ -M ./preliminarySubset_popmap -t 28 --fstats --vcf --structure --phylip --genepop --write-single-snp 
