#!/bin/bash

#----------------------------------------------------
# %%%%% PULLING IN DATA, CREATING SAMPLE SHEET %%%%%
#----------------------------------------------------
# Pull data from SNPsaurus sharesite (411 GB)
#wget -r -nH --cut-dirs=1 --no-parent --reject "index.html*" https://gc3fstorage.uoregon.edu/HG3YYDSX2/GC3F-PE-4637/

# Random stuff: extracting parts of the file names to create a sample sheet (use sed and awk for text matching)
#ls | cut -f1 -d"_" >> prelim_Barcodes.text
#sed -n '1~2!p' prelim_Barcodes.txt >> Barcodes.txt
#ls | awk -F'_' '{print $2}' >> prelim_SampleNames.txt
#sed -n '1~2!p' prelim_SampleNames.txt >> SampleNames.txt

#-----------------------------
# %%%%% FASTQC ANALYSIS %%%%%
#-----------------------------
# FastQC checks for the expected number of reads, read length, overall quality, and duplicate rate
# Using the CLI, FastQC can parse fastq.gz files. We'll get aggregate scores of all the sequences passed to the program
# Run this command in the directory containing the raw data (the .fastq.gz files)
# zcat on its own prints the output of unzipping these files to screen; the pipe directs the print output to fastqc, so files remain compressed
# Wildcards (?) are to prevent inclusion of sampples with failed barcodes (named "Undetermined.fastq.gz"); 8 is the number of cores specified

#zcat ????????-????????_*.fastq.gz | fastqc stdin:fastqc_results -t 8

#----------------------------------------
# %%%%% SETTING MAXIMUM FILE LIMITS %%%%%
#---------------------------------------
# For SNPsaurus NextRAD data, we have 382 pairs of sequences. Each of these pairs will generate 4 files in Stacks (Read 1, Read 2, rem 1, rem 2)
# Accessing that many files in a Linux system at once surpasses maximum open file limit, triggering an obscure error (From Stacks: "Error opening output file")
# Need to increase maximum allowable number of open files. Steps numbered below

# 1. First, check the limits on the number of open files
# "Hard" open file limit
#ulimit -Hn
# "Soft" open file limit
#ulimit -Sn

# 2. Two files need to be changed. Make copies of each file before altering originals (i.e. sudo cp originalFile backupFile)
# A. Access the /etc/security/limits.conf file, and add the relevant line

#sudo nano /etc/security/limits.conf
# Add below line to increase soft file limit
#user		soft	nofile	4096

# (Where 4096 is the amount you want to be allowable; increase to a maximum of the hard limit)
# (Tip: Lines in this file are of the format below)
#<domain/user>		<type>	<item>	<value>

# B. Access the /etc/pam.d/common-session file, then add the relevant line

# sudo nano /etc/pam.d/common-session
# Add below line
#session required pan_limits.so

# 3. Reboot computer

# 4. Recall ulimit commands to check limits successfully changed
#ulimit -Hn
#ulimit -Sn

#------------------------------
# %%%%% STACKS PROCESSING %%%%%
#-------------------------------
# Data is already demultiplexed. Use process_radtags to rename files, clean up messy sequences, and merge pair ends prior to alignment
# --c cleans the data, removing any read with an uncalled base; --q performs quality filtering: within a moving window (of length 0.15*total read length), if average score is Phred 10 or less, read is discarded
# --paired indicates paired end reads; -b specifies barcodes file (for renaming samples to original names); --index-index indicates that barcodes are in the headers of the FASTQ files (for forward and reverse reads)

# Preliminary runs showed the presence of Nextera transposase I7 (in R1 reads) and I5 (in R2 reads) sequences in NextRAD data.
# By looking at catalog of Nextera sequences from bbmap (bbmap/resources/nextera.fa.gz), the sequences for these adaptors were determined. These were passed to Stacks, for trimming.
# We specify a barcode mismatch of 2, in order to account for both the Nextera I5 sequence and the NextRAD I5 sequence (I7 sequence unchanged)

#process_radtags -p /RAID1/IMLS_GCCO/RawData/GC3F-PE-4637 --paired -b /RAID1/IMLS_GCCO/Analysis/Stacks/Barcodes.txt --index-index \
#-o /RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags -c -q --disable_rad_check --len-limit 100

process_radtags -p /RAID1/IMLS_GCCO/RawData/GC3F-PE-4637 --paired -b /RAID1/IMLS_GCCO/Analysis/Stacks/Barcodes.txt --index-index \
-o /RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags -c -q --adapter_1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --adapter_2 CTGTCTCTTATACACATCTGACGCTGCCGACGA --adapter_mm 2 --disable_rad_check --len-limit 100

# When process_radtags completes, separate the remainder reads into a new directory
# These reads represent samples that have fallen out of phase from the paired end read, due to a quality control removal
#mkdir RemainderReads
#mv *.rem.?.fq.gz RemainderReads/

#---------------------------------
# %%%%% REFERENCE ALIGNMENT %%%%%
#---------------------------------

# To align RAD loci with a reference, download and install GSNAP (Genomic Short-read Nucleotide Alignment Program, version 2021-07-23)
#wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-07-23.tar.gz

# Source Q. lobata reference genome from UCLA Valley Oak Genome Project page (and give it a discernible name) 
#wget -O Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz https://ucla.box.com/shared/static/9flhg1f4baygs7fb7hahz3knedlncvwd.gz
# Extract FASTA file (not necessary, as gmap_build has its own gunzip parameter)
#gunzip -k Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz

# Build Q. lobata reference index from FASTA file using gmap_build
#gmap_build -t 8 -d REF_Q.lobata_UCLA Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta

# Align output from process_radtags using gsnap. This command will interpret files as pairs (i.e. first two files will group together, then next two, then next two)
#16 cores; processing on compressed files; outputting SAM files

# Quercus acerifolia
#gsnap -t 16 --gunzip -d /ReferenceGenome/REF_Q.lobata_UCLA -A sam QUAC/*.fq.gz > QUAC_20210920.sam
# Quercus boyntonii
#gsnap -t 24 --gunzip -d /ReferenceGenome/REF_Q.lobata_UCLA -A sam QUBO/*.fq.gz > QUBO_20210922.sam
