#!/bin/bash

# This script runs the Stacks populations module, for each QUAC de novo assembly.
# It specifies the R80 argument, to filter missing dta (only loci shared between 80% of all samples are retained)

while IFS=, read -r assembly; do
	echo $assembly
	# Make an output directory for the R80 analysis
	mkdir ./output/$assembly/pop_R80
	# Call populations, using assembly variable to specify in path directory
        populations -P ./output/$assembly -O ./output/$assembly/pop_R80 -M ./popmap_QUAC_paramOpt -t 28 -R 80
done < ./output/QUAC_assemblies
