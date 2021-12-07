#!/bin/bash

# Script for extracting 4 metrics from Stacks output, in order to assess quality of denovo assembly for a specific parameter set

while IFS=, read -r assembly; do
	echo $assembly
        # Move into assembly folder
	cd ../output/$assembly

        # %%%%%%%%%%%%%%%%%%%%%%%
	# %%% EXTRACT METRICS %%%
	# %%%%%%%%%%%%%%%%%%%%%%%
	# COVERAGE
	stacks-dist-extract ./denovo_map.log cov_per_sample | cut -f 2 > metric-depth_of_cov
	# Consider: do we need to grab effective coverage too? Figure out after having conversation with Sean and Emily...

	# NUMBER OF ASSEMBLED LOCI
	stacks-dist-extract ./pop_R80/populations.log.distribs loci_per_sample | cut -f 2 > metric-assembled_loci

	# NUMBER OF POLYMORPHIC LOCI
	stacks-dist-extract ./pop_R80/populations.log.distribs snps_per_loc_postfilters | tail -n +4 | cut -f 2 > metric-polymorphic_loci

	# NUMBER OF SNPS
	stacks-dist-extract ./pop_R80/populations.log.distribs variant_sites_per_sample | cut -f 2 > metric-number_of_SNPs

	# PCR DUPLICATION RATE
	stacks-dist-extract ./gstacks.log.distribs effective_coverages_per_sample | cut -f 8 > metric-pcr_duplication_rate

	# %%%%%%%%%%%%%%%%%%%%
        # %%% BUILD MATRIX %%%
        # %%%%%%%%%%%%%%%%%%%%
	# Run R script to build a matrix in R containing metrics just extracted for the current assembly. The script will save this matrix to a .Rdata file
	Rscript /home/user/Documents/SSRvSNP/Code/denovoAssembly/parameterOptimization/build_assemblyMetricsMatrix.R

	# Move back into analysis folder
	cd ../../analysis

done < ../output/QUAC_assemblies
#done < ../output/QUAC_assemblies_subset
