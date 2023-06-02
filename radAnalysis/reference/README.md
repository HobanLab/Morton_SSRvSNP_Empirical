# Reference alignment

This folder contains the BASH scripts used to obtain and index reference genomes, align samples to reference genomes, 
and run the Stacks [ref_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) 
command used to generate loci datasets from reference alignments. 

Scripts are divided into 2 folders, one used for building reference alignments and one for SNP calling, as outlined below.

## Software
For indexing reference genomes and generating alignments, we used the [GMAP/GSNAP software](http://research-pub.gene.com/gmap/).
We analyzed reference alignment quality using the [samtools](http://www.htslib.org/download/) software.
We used the Stacks (v2.59) [ref_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) command for 
processing reference alignments, and the [populations module](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) for processing the SNP loci.

## Folders
### refAlignment
This folder contains the scripts we used for downloading and index reference genomes, and then building and processing reference alignments.
Note that the sample lists in the QUAC/QUBO folders contain more samples than those used for downstream analyses (i.e. replicates, Kessler individuals, etc.).

#### Reference genomes
For aligning _Quercus acerifolia_ (QUAC) samples, we used the [_Q. rubra_ reference genome](https://data.jgi.doe.gov/refine-download/phytozome?organism=Qrubra&expanded=687) (v2.0) available on Phytozome (preprint [here](https://www.authorea.com/users/618769/articles/643714-haplotype-resolved-chromosome-scale-genome-assembly-of-quercus-rubra-l?commit=52ec525eac731a3d3e1de15dd8d00b691deb4ee6)). 
For aligning _Q. boyntonii_ (QUBO) samples, we used the [_Q. robur_ reference genome](https://www.oakgenome.fr/?page_id=587) (V2_2N; formal publication [here](http://dx.doi.org/10.1038/s41477-018-0172-3)).

### processLoci
This folder contains the Stacks scripts and popmap files used to call SNPs from reference alignments and then process those SNP libraries (using the `populations` module).

For QUAC, we had to trim down the number of loci generated in our reference alignment datasets, in order to make them computationally feasible for STRUCTURE analyses.
We include the simple BASH script used to do this.
