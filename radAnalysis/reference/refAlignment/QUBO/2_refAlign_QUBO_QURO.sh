#!/bin/bash

# This script generates alignments between Quercus boyntonii NextRAD samples and the Q. robur reference genome using
# the GSNAP program ("promiscuous" alignment parameters from Paris et al. 2017, but with 4 allowed mismatches instead of 5).
# It also uses the samtools flagstat command to generate alignment statistics for each sample.

# Loop over the names in the list of QUBO samples (removing first row, column headers)
sed '1d' /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/QUBO_sampleList | cut -f1 |
while read sample; do
        # Create variables for each sample
        fq1=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.1.fq.gz # Forward reads
        fq2=/RAID1/IMLS_GCCO/Analysis/Stacks/process_RADtags/QUBO/${sample}.2.fq.gz # Reverse reads
        # Create sample name variables
        sam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP4/${sample}.sam
        bam=/RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP4/${sample}.bam
        # Align reads and process alignments. Indel penalty of 2, and prohibit terminal alignments through min-coverage=0.95 parameter (Paris et al. 2017), but fewer mismatches (4 instead of 5)
        gsnap -t 24 --gunzip -d ReferenceGenome/Q_robur/gsnap/REF_QURO_GMAP -A sam --omit-softclipped -i 2 -m 4 --min-coverage=0.95 $fq1 $fq2 > $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 24 $sam | samtools sort --threads 24 -o $bam
        # Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP4/QUBO_filtered_algnSummary_GSNAP4.tsv
        samtools flagstat -O tsv $bam >> /RAID1/IMLS_GCCO/Alignment/filteredReads/QUBO/GSNAP4/QUBO_filtered_algnSummary_GSNAP4.tsv
done
