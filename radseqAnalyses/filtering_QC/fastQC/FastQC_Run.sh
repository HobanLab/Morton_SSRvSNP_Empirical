#!/bin/bash

# FastQC checks for the expected number of reads, read length, overall quality, and duplicate rate
# Using the CLI, FastQC can parse fastq.gz files. We'll get aggregate scores of all the sequences passed to the program
# Run this command in the directory containing the raw data (the .fastq.gz files)
# zcat on its own prints the output of unzipping these files to screen; the pipe directs the print output to fastqc, so files remain compressed
# Wildcards (?) are to prevent inclusion of sampples with failed barcodes (named "Undetermined.fastq.gz"); 8 is the number of cores specified

# For processing raw, demultiplexed files (prior to running process_radtags)
zcat ????????-????????_*.fastq.gz | fastqc stdin:fastqc_results -t 8

# For trimmed (or untrimmed) labeled files, which have gone through the process_radtags step
zcat *.fq.gz | fastqc stdin -t 8
