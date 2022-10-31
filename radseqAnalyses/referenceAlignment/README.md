## REFERENCE ALIGNMENT

This folder contains the BASH scripts used to obtain and index reference genomes, align samples to reference genomes, 
and run the Stacks commands used to generate loci datasets from reference alignments. 

Scripts are divided into 2 folders, one for each species; the scripts and files in each folder are functionally identical,
and differ in the names of the directories (and in obtaining the two different reference genomes).

### Reference genomes
For aligning _Quercus acerifolia_ (QUAC) samples, we used the [_Q. rubra_ reference genome](https://data.jgi.doe.gov/refine-download/phytozome?organism=Qrubra&expanded=687) (v2.0) available on Phytozome. 
For aligning _Q. boyntonii_ (QUBO) samples, we used the [_Q. robur_ reference genome](https://www.oakgenome.fr/?page_id=587) (V2_2N).

### Software
For indexing reference genomes and generating alignments, we used the [GMAP/GSNAP software](http://research-pub.gene.com/gmap/).
We analyzed reference alignment quality using the [samtools](http://www.htslib.org/download/) software.
We used [Stacks](https://catchenlab.life.illinois.edu/stacks/) (v2.59) for processing reference alignments, and running the `populations` module.

### Outline
Within each folder, there are 3 BASH scripts and 2 population map files (popmaps). The popmap files are referenced in the scripts, and contain the samples analyzed 
in different steps. The two popmap files resemble each other, except the "_GardenWild" popmap files do not include duplicate samples.

#### Scripts
1. downloadAndIndex: these scripts contain the commands used for downloading the reference genomes and indexing those references (using the `gmap_build` command)
2. referenceAlignment: these scripts contain the commands used for aligning the NextRAD sequeqneces (after passing them through the Stacks `process_radtags` command) to their corresponding reference genomes.
These scripts also include the samtools flagstat command, for processing reference alignments (prior to analysis using Stacks)
3. _Stacks: these scripts contain the `ref_map.pl` and `populations` commands used from Stacks to generate loci datasets from reference alignments.
