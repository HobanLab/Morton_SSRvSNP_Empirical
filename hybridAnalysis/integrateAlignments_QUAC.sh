#!/bin/bash

# Script for running the stacks-integrate-alignments script for the BWA alignment between QUAC de novo consensus reads catalog and the Q. robur reference genome
# This integrates alignment positions into the de novo catalog, allowing more divergent loci to assembly while being able to provide information about loci colocalization
# Specify a specific temp directory. This is to avoid the error "No space left on device" when trying to write to a different tmp directory
export TMPDIR=/RAID1/TMP/

# Two filters are used: minimum alignment coverage (--min_alncov) and minimum percent identity (--min_pctid). Both are set to 80%
# These parameters and values were pulled from the Rivera-Colon and Catchen 2021 preprint
stacks-integrate-alignments --in-path /RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/ --bam-path /RAID1/IMLS_GCCO/Alignment/denovoConsensus/QUAC/BWA/QUAC_denovo_catalog_BWA.bam -O ./output/ --min_alncov 0.8 --min_pctid 0.8 --verbose

# After the new catalog has been constructed, run the populations program within the new directory, specifying the R80 parameter (as in de novo assemblies)
populations -P ./output -O ./populations --popmap ./QUAC_popmap --threads 28 -R 80 
