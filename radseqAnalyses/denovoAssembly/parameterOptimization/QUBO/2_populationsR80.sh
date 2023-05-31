#!/bin/bash

# This script runs the Stacks populations module, for each QUBO de novo assembly.
# It specifies the R80 argument, to filter missing dta (only loci shared between 80% of all samples are retained)

while IFS=, read -r assembly; do
	echo $assembly
	# Make an output directory for the updated analysis
	mkdir ./output/$assembly/updated
	# Call populations, using assembly variable to specify in path directory
        populations -P ./output/$assembly -O ./output/$assembly/updated -M ./popmap_QUBO_paramOpt -t 28 -R 80
	sleep 5s
done < ./output/QUBO_assemblies
