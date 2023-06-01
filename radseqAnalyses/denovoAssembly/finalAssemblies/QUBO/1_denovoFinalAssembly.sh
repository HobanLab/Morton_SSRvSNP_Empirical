#!/bin/bash

# Script for constructing the finalized de novo assembly for Quercus boyntonii samples
# m, M, n, and gt-alpha were determined through parameter optimization; N of 0 is suggested in Paris et al. 2019 for assemblies with high coverage; R -80 is to compare samples present in 80% of individuals
denovo_map.pl -T 28 -m 7 -M 5 -n 5 --gt-alpha 0.01 -N 0 --samples /RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO --popmap ./QUBO_popmap -o ./output/ --paired --rm-pcr-duplicates
