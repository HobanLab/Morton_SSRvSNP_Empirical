#!/bin/bash

# This script contains the Stacks ref_map.pl command for filtered, reference aligned QUBO reads.
# It also analyzes the outputs of the ref_map step, using the stacks-dist-extract tool

# %%% IDENTIFYING REFERENCE ALIGNMENT LOCI %%%
# The only modules that run for ref_map.pl are gstacks and populations;
# we pass the --rm-pcr-duplicates and --gt-alpha arguments to gstacks in order to mirror the parameters used for de novo assembly.
# Paired end reads were passed to the GMAP software, so Stacks should automatically detect that it's working with paired end data

# Command for processing reference alignments, specifying the input (--samples) and output directories
ref_map.pl -T 28 --gt-alpha 0.01 --rm-pcr-duplicates --samples /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP4 --popmap ./QUBO_popmap -o ./output/

# After ref_map.pl script has completed, analyze the gstacks outputs to assess reference alignment quality
stacks-dist-extract gstacks.log.distribs bam_stats_per_sample >> bamStats.csv
