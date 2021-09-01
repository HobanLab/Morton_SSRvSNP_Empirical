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
# Run this command in the directory containing the raw data (i.e. .fastq.gz files). 
# zcat on its own prints the output of unzipping these files to screen; the stdin portion takes that output and pushes it to FastQC
# test_results is the name of the folder; 8 is the number of cores specified

#zcat *fastq.gz | fastqc stdin:test_results -t 16

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

# 2. Two files need to changed. Make copies of each file before altering originals (i.e. sudo cp originalFile backupFile)
# FIRST FILE--access the file, then add the relevant line

#sudo nano /etc/security/limits.conf
# Add below line to increase soft file limit
#user		soft	nofile	4096

# (Where 4096 is the amount you want to be allowable; increase to a maximum of the hard limit)
# (Tip: Lines in this file are of the format below)
#<domain/user>		<type>	<item>	<value>

# 3. SECOND FILE--access the file, then add the relevant line

# sudo nano /etc/pam.d/common-session
# Add below line
#session required pan_limits.so

# 4. Reboot computer

# 5. Recall ulimit commands to check limits successfully changed
#ulimit -Hn
#ulimit -Hn

#------------------------------
# %%%%% STACKS PROCESSING %%%%%
#-------------------------------
# Data is already demultiplexed. Use process_radtags to clean up messy sequences, and to merge pair ends.
# --c cleans the data, removing any read with an uncalled base. --q performs quality filtering: within a moving window (of length 0.15*total read length), if average score is Phred 10 or less, read is discarded
# --paired indicates paired end reads
#process_radtags -p /RAID1/IMLS_GCCO/RawData/GC3F-PE-4637 --paired -o /RAID1/IMLS_GCCO/Analysis/Stacks/ProcessRADtags -c -q --disable_rad_check --len-limit 100

# New command, with barcodes to rename files
process_radtags -p /RAID1/IMLS_GCCO/RawData/GC3F-PE-4637 --paired -b /RAID1/IMLS_GCCO/Analysis/Stacks/Barcodes.txt --index-index \
-o /RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags -c -q --disable_rad_check --len-limit 100

#---------------------------------
# %%%%% REFERENCE ALIGNMENT %%%%%
#---------------------------------

# To align RAD loci with a reference, using GSNAP (Genomic Short-read Nucleotide Alignment Program, version 2021-07-23)
#wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-07-23.tar.gz

# Source Q. lobata reference genome from UCLA Valley Oak Genome Project page (and give it a discernible name) 
#wget -O Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz https://ucla.box.com/shared/static/9flhg1f4baygs7fb7hahz3knedlncvwd.gz
# Extract FASTA file (not necessary, as gmap_build has its own gunzip parameter
#gunzip -k Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz

# Build Q. lobata reference index from FASTA file using gmap_build
#gmap_build -t 8 -d REF_Q.lobata_UCLA Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta

# Align output from process_radtags using gsnap. This command will interpret files as pairs (i.e. first two files will group together, then next two, then next two)
#gsnap --gunzip -t 16 -D /RAID1/IMLS_GCCO/ReferenceGenome/REF_Q.lobata_UCLA <ALL_FASTA_FILES_HERE>

