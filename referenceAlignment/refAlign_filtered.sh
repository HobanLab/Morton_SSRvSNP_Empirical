#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%% REFERENCE ALIGNMENT -- FILTERED READS %%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script includes commands used to obtain reference genomes and align reads (filtered using the Stacks process_radtags command) with those reference genomes, for the IMLS_GCCO RADseq vs. microsatellite project

# %%%%%%%%%%%%%%%%%%%%
# %%%%% SOFTWARE %%%%%
# %%%%%%%%%%%%%%%%%%%%

# The BWA (Burrows-Wheeler Aligner) is an open source software that allows for short, deep-sequencing reads to be aligned to long reference sequences.
# It is available at https://sourceforge.net/projects/bio-bwa/files/

# GSNAP is a tool to align single- and paired-end reads to a reference genome, and works for sequences as short as 14 nt. Unlike BWA, GSNAP allows users to prohibit soft-clipping, which is why it was utilized for Stacks processes.
# It is available at http://research-pub.gene.com/gmap/

# samtools is a software package comprised of tools for manipulating alignment (.sam) and compressed alignment (.bam) files. The flagstat command from samtools is used in this script to summarize alignment qualities.
# It is available at http://www.htslib.org/download/

# Commands used to process these alignments (using the Stacks ref_map.pl script) are in other scripts

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% OBTAIN AND INDEX QUERCUS ROBUR REFERENCE GENOME %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%% BWA %%%%
# Source Q. robur reference genome assembly (V2_2N) from French consortium led by Christophe Plomion at INRA Bordeaux (and give it a discernible name)
wget -O Qrobur_V2_2N.fa.gz https://urgi.versailles.inra.fr/download/oak/Qrob_V2_2N.fa.gz
# Unzip
gunzip Qrobur_V2_2N.fa.gz
# Create index using BWA. This index is what allows genes to be mapped to chromosome regions at the end of an analysis
bwa index Qrobur_V2_2N.fa -p Qrobur_db

# %%%% GSNAP %%%%
# Build Q. robur reference index from FASTA file using gmap_build
gmap_build -t 8 -d REF_QURO_GSNAP Qrobur_V2_2N.fa

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% CREATE ALIGNMENTS %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%% QUAC %%%%
# The blocks below contain similar commands, but are iterated for the 3 different alignment programs/parameters explored: BWA (default), GSNAP (prohibit soft-clipping), and GSNAP2 ("promiscuous" alignment parameters)

# ---BWA---
# Loop over the sample names listed in the popmap file for QUAC samples (removing first row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/QUAC_popmap | cut -f1 |
while read sample; do
	#Create variables for each sample
        fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.1.fq.gz # Forward reads
        fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.2.fq.gz # Reverse reads
        # Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/${sample}.bam
        # Align reads and process alignments
        bwa mem /RAID1/IMLS_GCCO/ReferenceGenome/Q_robur/Qrobur_db $fq1 $fq2 -t 24 -o $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
	# Generate basic alignment statistics using flagstat, and output them to a summary file
	echo ${sample} >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/QUAC_filtered_algnSummary_BWA.tsv
	samtools flagstat -O tsv $bam >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/QUAC_filtered_algnSummary_BWA.tsv
done

# ---GSNAP---
# Loop over the sample names listed in the popmap file for QUAC samples (removing first row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/QUAC_popmap | cut -f1 |
while read sample; do
        #Create variables for each sample
        fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.1.fq.gz # Forward reads
        fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.2.fq.gz # Reverse reads
        # Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP/${sample}.bam
        # Align reads and process alignments
	gsnap -t 24 --gunzip -d ReferenceGenome/Q_robur/gsnap/REF_QURO_GMAP -A sam --omit-softclipped $fq1 $fq2 > $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
        # Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP/QUAC_filtered_algnSummary_GSNAP.tsv
        samtools flagstat -O tsv $bam >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP/QUAC_filtered_algnSummary_GSNAP.tsv
done

# ---GSNAP2---
# Loop over the sample names listed in the popmap file for QUAC samples (removing first row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/QUAC_popmap | cut -f1 |
while read sample; do
        #Create variables for each sample
        fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.1.fq.gz # Forward reads
        fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUAC/${sample}.2.fq.gz # Reverse reads
        # Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP2/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP2/${sample}.bam
        # Align reads and process alignments. Allow 5 mismatches ("promiscuous"), indel penalty of 2, and prohibit terminal alignments through min-coverage=0.95 parameter (Paris et al. 2017)
        gsnap -t 16 --gunzip -d ReferenceGenome/Q_robur/gsnap/REF_QURO_GMAP -A sam --omit-softclipped -m 5 -i 2 --min-coverage=0.95 $fq1 $fq2 > $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
        # Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP2/QUAC_filtered_algnSummary_GSNAP2.tsv
        samtools flagstat -O tsv $bam >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUAC/GSNAP2/QUAC_filtered_algnSummary_GSNAP2.tsv
done

# %%%% QUBO %%%%
# The blocks below contain similar commands, but are iterated for the 2 different alignment programs/parameters explored: BWA (default) and GSNAP (prohibitting soft-clipping)

# ---BWA---
# Loop over the sample names listed in the popmap file for QUBO samples (removing firts row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/QUBO_popmap | cut -f1 |
while read sample; do
        #Create variables for each sample
        fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.1.fq.gz # Forward reads
        fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.2.fq.gz # Reverse reads
	# Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/${sample}.bam
        # Align reads and process alignments
        bwa mem /RAID1/IMLS_GCCO/ReferenceGenome/Q_robur/Qrobur_db $fq1 $fq2 -t 24 -o $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
	# Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/QUBO_filtered_algnSummary_BWA.tsv
        samtools flagstat -O tsv $bam >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/QUBO_filtered_algnSummary_BWA.tsv
done

# ---GSNAP---
# Loop over the sample names listed in the popmap file for QUBO samples (removing first row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/QUBO_popmap | cut -f1 |
while read sample; do
        #Create variables for each sample
        fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.1.fq.gz # Forward reads
        fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.2.fq.gz # Reverse reads
        # Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP/${sample}.bam
        # Align reads and process alignments
        gsnap -t 24 --gunzip -d ReferenceGenome/Q_robur/gsnap/REF_QURO_GMAP -A sam --omit-softclipped $fq1 $fq2 > $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
        # Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP/QUBO_filtered_algnSummary_GSNAP.tsv
        samtools flagstat -O tsv $bam >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP/QUBO_filtered_algnSummary_GSNAP.tsv
done

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% OUTDATED: OBTAIN AND INDEX QUERCUS LOBATA REFERENCE GENOME %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Source Q. lobata reference genome from UCLA Valley Oak Genome Project page (and give it a discernible name) 
wget -O Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz https://ucla.box.com/shared/static/9flhg1f4baygs7fb7hahz3knedlncvwd.gz
# Extract FASTA file (not necessary, as gmap_build has its own gunzip parameter)
gunzip -k Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz

# Build Q. lobata reference index from FASTA file using gmap_build
gmap_build -t 8 -d REF_Q.lobata_UCLA Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta

# Align output from process_radtags using gsnap. This command will interpret files as pairs (i.e. first two files will group together, then next two, then next two)
#16 cores; processing on compressed files; outputting SAM files

# Quercus acerifolia
gsnap -t 24 --gunzip -d /ReferenceGenome/REF_Q.lobata_UCLA -A sam QUAC/*.fq.gz > QUAC_20210928.sam
# Quercus boyntonii
gsnap -t 24 --gunzip -d /ReferenceGenome/REF_Q.lobata_UCLA -A sam QUBO/*.fq.gz > QUBO_20210929.sam
