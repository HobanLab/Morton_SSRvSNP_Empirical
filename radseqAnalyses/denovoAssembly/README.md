## _De novo_ Assemblies

This folder contains the scripts and files used to run _de novo_ assembly parameter optimization, to analyze the outputs
of the parameter optimizaiton, and to run the [denovo_map.pl command](https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) 
for the optimized _de novo_ assemblies. The contents of the two folders are outlined below. 

##### denovo.params
This text file is referenced for both species, and is simply a list of all the possible combinations for the parameters varied
in the _de novo_ assembly.

#### plot_deNovoAssemblyMetrics_MS.R
This script generates plots that summarize the assembly metrics used to analyze parameter optimizations (i.e.
number of loci/SNPs, number of polymorphic loci, coverage, etc.). The plots in this script were included in the supplemental
material for the manuscript supplement

### parameterOptimization
This folder contains the scripts used for _de novo_ assembly parameter optimization. It is split into two folders, one for each
species; the contents for each folder are largely identical, but with the appropriate name for each species, and are numbered 
according to the order that they're called in the parameter optimization process. 

### finalAssemblies
This folder contains the [denovo_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) calls for the finalized
_de novo_ assemblies, for each species. 

