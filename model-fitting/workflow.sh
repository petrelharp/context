#!/bin/bash

# simulate up some sequence for testing
Rscript ../sim-seq.R -c cpg-model.json -t .1 -s 1000 -d simseqs -o sim-cpg-123456.RData
# and count the paired tuples, width 4 = 1+2+1
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-4 with all one-base transitions
Rscript ../flex-genmat.R -c genmatrices/singlebase.json -w 4
#   width-5 with all one-base transitions
Rscript ../flex-genmat.R -c genmatrices/singlebase.json -w 5
#   width-4 with all two-base transitions
Rscript ../flex-genmat.R -c genmatrices/dualbases.json -w 4

# fit a simple model
Rscript ../fit-model.R -i sim-cpg-123456.4.2.l1.counts -l 1 -m genmatrices/genmatrix-4-singlebase.RData -o simple-fit.RData
# and a more complex one
Rscript ../fit-model.R -i sim-cpg-123456.4.2.l1.counts -l 1 -m genmatrices/genmatrix-4-dualbases.RData -o complex-fit.RData

# look at residuals, width 4
Rscript ../compute-resids.R -i sim-cpg-123456-modelfit-54321.RData

# look at residuals, width 5
#   first compute counts
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 5 -s 3 -l 1
Rscript ../compute-resids.R -i sim-cpg-123456-modelfit-54321.RData
