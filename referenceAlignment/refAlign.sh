#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%% REFERENCE ALIGNMENT %%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script includes commands used to obtain reference genomes and align samples with those reference genomes, for the IMLS_GCCO RADseq vs. microsatellite project
# Many of these steps utilize the BWA (Burrows-Wheeler Aligner), which is an open source software that allows for short, deep-sequencing reads to be aligned to long reference sequences.
# It is available at https://sourceforge.net/projects/bio-bwa/files/
# Commands used to process these alignments (using Stacks) are in other scripts

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% OBTAIN QUERCUS ROBUR REFERENCE GENOME %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Source Q. robur reference genome assembly (V2_2N) from French consortium led by Christophe Plomion at INRA Bordeaux (and give it a discernible name)
#wget -O Qrobur_V2_2N.fa.gz https://urgi.versailles.inra.fr/download/oak/Qrob_V2_2N.fa.gz
# Unzip
#gunzip Qrobur_V2_2N.fa.gz
# Create index using BWA. This index is what allows genes to be mapped to chromosome regions at the end of an analysis
#bwa index Qrobur_V2_2N.fa -p Qrobur_db

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% CREATE ALIGNMENTS %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%% QUAC %%%%
# Loop over the sample names listed in the popmap file for QUAC samples (removing first row, column headers)
#sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/QUAC_popmap | cut -f1 |
#while read sample; do
	#Create variables for each sample
        #fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.1.fq.gz # Forward reads
        #fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.2.fq.gz # Reverse reads
        # Create sample name variables
        #sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/${sample}.sam
        #bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/${sample}.bam
        # Align reads and process alignments
        #bwa mem /RAID1/IMLS_GCCO/ReferenceGenome/Q_robur/Qrobur_db $fq1 $fq2 -t 24 -o $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        #samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
	# Generate basic alignment statistics using flagstat, and output them to a summary file
	#echo ${sample} >> QUAC_filtered_algnSummary.txt
	#samtools flagstat $bam >> QUAC_filtered_algnSummary.txt
#done

# %%%% QUBO %%%%
# Loop over the sample names listed in the popmap file for QUBO samples (removing firts row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/QUBO_popmap | cut -f1 |
while read sample; do
        #Create variables for each sample
        #fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.1.fq.gz # Forward reads
        #fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.2.fq.gz # Reverse reads
	# Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/${sample}.bam
        # Align reads and process alignments
        #bwa mem /RAID1/IMLS_GCCO/ReferenceGenome/Q_robur/Qrobur_db $fq1 $fq2 -t 24 -o $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        #samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
	# Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> QUBO_filtered_algnSummary.txt
        samtools flagstat $bam >> QUBO_filtered_algnSummary.txt
done

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% OUTDATED: OBTAIN QUERCUS LOBATA REFERENCE GENOME %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Source Q. lobata reference genome from UCLA Valley Oak Genome Project page (and give it a discernible name) 
#wget -O Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz https://ucla.box.com/shared/static/9flhg1f4baygs7fb7hahz3knedlncvwd.gz
# Extract FASTA file (not necessary, as gmap_build has its own gunzip parameter)
#gunzip -k Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz

# Build Q. lobata reference index from FASTA file using gmap_build
#gmap_build -t 8 -d REF_Q.lobata_UCLA Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta

# Align output from process_radtags using gsnap. This command will interpret files as pairs (i.e. first two files will group together, then next two, then next two)
#16 cores; processing on compressed files; outputting SAM files

# Quercus acerifolia
#gsnap -t 24 --gunzip -d /ReferenceGenome/REF_Q.lobata_UCLA -A sam QUAC/*.fq.gz > QUAC_20210928.sam
# Quercus boyntonii
#gsnap -t 24 --gunzip -d /ReferenceGenome/REF_Q.lobata_UCLA -A sam QUBO/*.fq.gz > QUBO_20210929.sam
