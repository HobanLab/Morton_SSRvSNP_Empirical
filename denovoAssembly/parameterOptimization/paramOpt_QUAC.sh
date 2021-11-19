#!/bin/bash

# Script for looping over parameters file, detailing the parameter space for the Stacks de novo assembly
# 4 variables (coming after the -r value) change values
while IFS=, read -r m M n gtalpha; do
	# Make output directory to store results, based on combination of parameter values
	mkdir ./output/QUAC_m${m}_M${M}_n${n}_ga${gtalpha}
	# m values range from 3--7; M/n range from 1 to 8; N is default (M+2) but, depending on optimized m value, will probably be set to 0; gt-alpha is 0.05 (more liberal) or 0.01 (stricter)
	echo 'Calling Stacks denovo_map using the following values--m:' "$m" ', M:'  "$M" ', n:'  "$n" ', gt-alpha:' "$gtalpha" &
	denovo_map.pl -T 28 -m $m -M $M -n $n --gt-alpha $gtalpha --samples ../../process_RADtags/QUAC/ --popmap ./popmap_QUAC_paramOpt -o ./output/QUAC_m${m}_M${M}_n${n}_ga${gtalpha}/ --paired --rm-pcr-duplicates
	sleep 5s
done < ../denovo.params
echo "All QUAC parameter optimization runs complete!"
