#!/bin/bash

# Preliminary analysis on a subset of samples, to ensure some basic expectations of our population genetic data
# -m is left to default (3); -M and -n are set equal to one another, but value is more or less arbitrary. 
# Remove PCR duplicates, generate Fst stats,  and multithread
denovo_map.pl -T 28 -M 4 -n 4 --samples ./ --popmap ./preliminarySubset_popmap -o ./output/ --paired --rm-pcr-duplicates
