#!/bin/bash

python ../tuple-counting/get-regions-from-alignment.py --axtfile mm9rn5/chrX.mm9.rn5.net.axt.gz --posfile <(tail -n +2 10_filter_regions.pos | cut -f 1 | sed -e 's/\(.*\):\([0-9]*\)-\([0-9]*\).*/\1 \2 \3/') -o mm9rn5/chrX.mm9.rn5.net.axt.sub.gz

python ../tuple-counting/count-paired-tuples.py --strict --infile mm9rn5/chrX.mm9.rn5.net.axt.sub.gz --informat axt -w 5 -l 1 -r 1 --outfile mm9rn5/chrX.mm9.rn5.5.1.counts
