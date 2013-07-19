#!/bin/bash

for x in $(ls -d -1 /home/mkessler/Hsap_ANALYSIS/23_peter_files/*);
do
    outfile=$(echo $x | sed -e "s/\.out\./.counts./" | sed -e "s/.*\///")
    # echo $outfile
    python /home/peter/projects/codons/count-tuples.py $x | gzip -c > $outfile
done

