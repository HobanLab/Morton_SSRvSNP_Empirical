#!/bin/bash

# This includes command used for analysis of samples aligned to a reference

#---------------------------------
# %%%%% REFERENCE ALIGNMENT %%%%%
#---------------------------------

# To align RAD loci with a reference, download and install GSNAP (Genomic Short-read Nucleotide Alignment Program, version 2021-07-23)
#wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-07-23.tar.gz

# %%%%%%%%%%%%%%%%%%%%%%%
# %%%% QUERCUS ROBUR %%%%
# %%%%%%%%%%%%%%%%%%%%%%%
# Source Q. robur reference genome assembly (V2_2N) from French consortium led by Christophe Plomion at INRA Bordeaux (and give it a discernible name)
#wget -O Qrobur_V2_2N.fa.gz https://urgi.versailles.inra.fr/download/oak/Qrob_V2_2N.fa.gz

# %%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS LOBATA %%%%
# %%%%%%%%%%%%%%%%%%%%%%%
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
