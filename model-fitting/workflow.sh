#!/bin/bash

# simulate up some sequence for testing
Rscript ../sim-seq.R -c cpg-model.json -t .01 -s 100 -d simseqs
# and count the paired tuples
Rscript ../count-seq.R -i simseqs/simseq-2014-08-01-14-18-0169113.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-5 with all one-base transitions
Rscript ../flex-genmat.R -c genmatrices/singlebase.mut -w 4
#   width-5 with all two-base transitions
Rscript ../flex-genmat.R -c genmatrices/dualbase.mut -w 4

# fit a simple model
