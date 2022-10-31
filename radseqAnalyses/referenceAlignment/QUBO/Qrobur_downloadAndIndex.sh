#!/bin/bash

# Commands used for downloading the Q. robur genome, and indexing that genome using the gmap_build command

# Source Q. robur reference genome assembly (V2_2N) from French consortium led by Christophe Plomion at INRA Bordeaux (and give it a discernible name)
wget -O Qrobur_V2_2N.fa.gz https://urgi.versailles.inra.fr/download/oak/Qrob_V2_2N.fa.gz

# Command for indexing the genome, which will create a genome "database" from a FASTA file. Specify 16 threads
gmap_build -t 16 -d REF_QURO_GSNAP Qrobur_V2_2N.fa
