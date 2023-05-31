#!/bin/bash

# FastQC checks for the expected number of reads, read length, overall quality, and duplicate rate
# Using the CLI, FastQC can parse fastq.gz files. We'll get aggregate scores of all the sequences passed to the program

# We ran this command in the directory containing the data (.fastq.gz or .fq.qz files).
# We analyzed 3 datasets with FasQC:
#   + "Prefiltering": raw data, before processing with process_radtags
#   + "Untrimmed": after processing with process_radtags, with no adapter trimming specified in process_radtags
#   + "Filtering": after processing with process_radtags, with adapter trimming specified in process_radtags
# Our de novo assemblies and reference alignments were built using the outputs from "Filtering"

# zcat on its own prints the output of unzipping these files to screen;
# the pipe directs the print output to fastqc, so files remain compressed
# Wildcards (?) are to prevent inclusion of sampples with failed barcodes (named "Undetermined.fastq.gz")
# 8 cores are specified.

# For processing raw, demultiplexed files (prior to running process_radtags; those in Prefiltering)
zcat ????????-????????_*.fastq.gz | fastqc stdin:fastqc_results -t 8

# For Trimmed (or Untrimmed) labeled files, which have gone through the process_radtags step
zcat *.fq.gz | fastqc stdin -t 8
