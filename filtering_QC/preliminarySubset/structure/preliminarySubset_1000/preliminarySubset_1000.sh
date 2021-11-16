#!/bin/bash/

for k in {2..7} ;
do
    for r in {1..50} ;
    do
        /home/user/Documents/SSRvSNP/software/structure/console/structure -i populations.structure -m mainparams.txt -e extraparams.txt -K $k -o output/Demo_output_$k-$r > output/outfile_$k_$r  &
        sleep 3s
    done
done
echo "All runs started"
