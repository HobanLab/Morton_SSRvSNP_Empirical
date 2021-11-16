#!/bin/bash

# The below commands are used to extract data from compressed (*.gz) files and pass it to fastqc, which will generate a report of data quality
#Each of these commands are sent within a directory containing the sample files.

# For processing raw, demultiplexed files (prior to running process_radtags)
zcat ????????-????????_*.fastq.gz | fastqc stdin:fastqc_results -t 8

# For trimmed (or untrimmed) labeled files, which have gone through the process_radtags step
zcat *.fq.gz | fastqc stdin -t 8
