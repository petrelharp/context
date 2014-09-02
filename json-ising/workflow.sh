#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c ising-model.json -t .1 -s 100000 -o test-ising-123456.RData

# and count the Tmers
Rscript ../count-seq.R -i test-ising-123456.RData -w 3 -s 1 -l 1

# precompute generator matrices:
#   width-3
Rscript ../make-genmat.R -c genmatrices/complete.json -w 3

# fit a model
Rscript ../fit-model.R -i test-ising-123456.3.1.l1.counts -l 1 -m genmatrices/genmatrix-3-complete.RData -j 54321

# mcmc also
Rscript ../mcmc-model.R -i test-ising-123456-genmatrix-3-complete-54321.RData -c ising-model.json -b 3 -j 1111
Rscript ../mcmc-model.R -i test-ising-123456-genmatrix-3-complete-54321-mcmc-1111.RData -b 100 -j 2222
Rscript ../mcmc-model.R -i test-ising-123456-genmatrix-3-complete-54321-mcmc-2222.RData -b 1000 -j 3333


# OK, now do this on more Tmer sizes:
Rscript ../count-seq.R -i test-ising-123456.RData -w 4 -s 2 -l 1
Rscript ../make-genmat.R -c genmatrices/complete.json -w 4
Rscript ../fit-model.R -i test-ising-123456.4.2.l1.counts -l 1 -m genmatrices/genmatrix-4-complete.RData -j 54321

Rscript ../count-seq.R -i test-ising-123456.RData -w 5 -s 3 -l 1
Rscript ../make-genmat.R -c genmatrices/complete.json -w 5
Rscript ../fit-model.R -i test-ising-123456.5.3.l1.counts -l 1 -m genmatrices/genmatrix-5-complete.RData -j 54321

Rscript ../count-seq.R -i test-ising-123456.RData -w 6 -s 2 -l 2
Rscript ../make-genmat.R -c genmatrices/complete.json -w 6
Rscript ../fit-model.R -i test-ising-123456.6.2.l2.counts -l 2 -m genmatrices/genmatrix-6-complete.RData -j 54321


###
# and for the constrained model, with longer windows
Rscript ../sim-seq.R -c ising-model-constrained.json -t .1 -s 1000000 -o constrained-ising-123.RData
Rscript ../count-seq.R -i constrained-ising-123.RData -w 9 -s 3 -l 3 -o constrained-ising-123.9.3.l3.counts
Rscript ../make-genmat.R -c ising-model-constrained.json -w 9 -o genmatrices/genmatrix-9-ising-model-constrained.RData
Rscript ../fit-model.R -i constrained-ising-123.9.3.l3.counts -l 3 -m genmatrices/genmatrix-9-complete.RData -o constrained-ising-123-genmatrix-9-ising-model-constrained.RData
Rscript ../mcmc-model.R -i constrained-ising-123-genmatrix-9-ising-model-constrained.RData -c ising-model-constrained.json -b 3 -o constrained-ising-123-genmatrix-9-mcmc-1.RData
Rscript ../mcmc-model.R -i constrained-ising-123-genmatrix-9-ising-model-constrained.RData -c ising-model-constrained.json -b 100 -o constrained-ising-123-genmatrix-9-mcmc-2.RData
