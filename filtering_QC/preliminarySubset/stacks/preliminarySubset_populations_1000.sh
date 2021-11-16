#!/bin/bash

# %%%%% PRELIMINARY ANALYSIS %%%%%
# We assembled loci de novo using a representative subset of all of our samples (both species, both garden and wild samples). This was to ensure some basic checks about our data:
# 1. Replicates should match one another
# 2. Samples from the same populations should resemble one another more than they do samples from other populations (for wild samples)

# Generate a whitelist of 1000 randomly selected loci, to analyze preliminary subset. This is to speed up analysis.
grep -v "^#" populations.sumstats.tsv | cut -f 1 | sort | uniq | shuf | head -n 1000 | sort -n > wl_1000.tsv

# Run standalone Stacks populations program, and generate VCF, STRUCTURE, genepop, and Fst outputs. Use the whitelist of 1000 loci and only analyze first SNP per locus
populations -P ./output/ -O ./output/populations/ -M ./preliminarySubset_popmap -t 28 --fstats --vcf --structure --phylip --genepop --write-single-snp -W ./output/wl_1000.tsv
