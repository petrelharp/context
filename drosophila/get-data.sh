#!/bin/bash

# Pairwise alignments: ( http://hgdownload.soe.ucsc.edu/goldenPath/dm3/vsDroSim1/ )
#  axtNet/*.dm3.droSim1.net.axt.gz: chained and netted alignments,
#  i.e. the best chains in the D. melanogaster genome, with gaps in the best
#  chains filled in by next-best chains where possible.  The axt format is
#  described in http://genome.ucsc.edu/goldenPath/help/axt.html

for chrom in 2L 2R 3L 3R 4 M X;
do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/vsDroSim1/axtNet/chr${chrom}.dm3.droSim1.net.axt.gz -P data
    wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/vsDroSec1/axtNet/chr${chrom}.dm3.droSec1.net.axt.gz -P data
done


