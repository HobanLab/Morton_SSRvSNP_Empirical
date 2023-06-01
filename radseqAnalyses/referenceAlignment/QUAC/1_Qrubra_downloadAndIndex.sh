#!/bin/bash

# Commands used for downloading the Q. rubra genome, and indexing that genome using the gmap_build command

# curl command used to download QURA reference genome (v2.0) from JGI website (Phytozome), along with README and Data Release Policy
curl --cookie jgi_session=/api/sessions/0e83085ab6485fb58e7642c9ac876443 --output download.20220919.094251.zip -d "{\"ids\":{\"Phytozome-687\":[\"60bfb62fc399d4ad32fd480d\",\"60bfb630c399d4ad32fd4811\",\"60bfb62fc399d4ad32fd4810\"]}}" -H "Content-Type: application/json" https://files.jgi.doe.gov/filedownload/

# Command for indexing the genome, which will create a genome "database" from a FASTA file. Specify 16 threads
gmap_build -D /RAID1/IMLS_GCCO/ReferenceGenome/Q_rubra -t 16 -d REF_QURU_GMAP Phytozome/PhytozomeV13/Qrubra/v2.1/assembly/Qrubra_687_v2.0.fa

