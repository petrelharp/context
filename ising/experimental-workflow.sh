#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../scripts/sim-seq.R -c ising-model.json -t .1 -s 100000 -d simseqs -o test-ising-123456.RData

# and count the Tmers
Rscript ../scripts/count-seq.R -i simseqs/test-ising-123456.RData -w 3 -s 1 -l 1
Rscript ../scripts/count-seq.R -i simseqs/test-ising-123456.RData -w 9 -s 5 -l 2

# precompute generator matrices:
#   width-3
Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 3
Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 4
Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 9

echo "check simulated model matches expected"
../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd simseqs/test-ising-123456.RData genmatrices/genmatrix-3-complete.RData 

echo "fitting the model"
Rscript ../scripts/fit-model.R -i simseqs/test-ising-123456-3-root-1-tip-l1.counts -m genmatrices/genmatrix-3-complete.RData -j 54321

echo "computing residuals"
Rscript ../scripts/compute-resids.R -i simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321.RData

echo "look at results"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321.RData simseqs/test-ising-123456.RData

 
echo "mcmc also"
Rscript ../scripts/mcmc-model.R -i simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321.RData -c ising-model.json -b 3 -j 1111
Rscript ../scripts/mcmc-model.R -i simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321.RData -c ising-model.json -b 100 -j 2222
Rscript ../scripts/mcmc-model.R -i simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321.RData -c ising-model.json -b 1000 -j 3333

echo "look at results"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321-mcmc-2222.RData simseqs/test-ising-123456.RData
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/test-ising-123456-3-root-1-tip-l1-genmatrix-3-complete-54321-mcmc-3333.RData simseqs/test-ising-123456.RData

echo "mcmc on offset counts"
Rscript ../scripts/count-seq.R -i simseqs/test-ising-123456.RData -w 4 -s 2 -l 1 --shift 4
Rscript ../scripts/fit-model.R -i simseqs/test-ising-123456-4-root-2-tip-l1-shift4.counts.1 -m genmatrices/genmatrix-4-complete.RData -j 10101
Rscript ../scripts/mcmc-model.R -i simseqs/test-ising-123456-4-root-2-tip-l1-shift4-genmatrix-4-complete-10101.RData -c ising-model.json -b 1000 -j 20202
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/test-ising-123456-4-root-2-tip-l1-shift4-genmatrix-4-complete-10101-mcmc-20202.RData simseqs/test-ising-123456.RData

echo "now on a larger T-mer (4/2)"
Rscript ../scripts/count-seq.R -i simseqs/test-ising-123456.RData -w 4 -s 2 -l 1
Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 4
Rscript ../scripts/fit-model.R -i simseqs/test-ising-123456-4-root-2-tip-l1.counts -m genmatrices/genmatrix-4-complete.RData -j 001
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/test-ising-123456-4-root-2-tip-l1-genmatrix-4-complete-001.RData simseqs/test-ising-123456.RData

# Rscript ../scripts/count-seq.R -i test-ising-123456.RData -w 5 -s 3 -l 1
# Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 5
# Rscript ../scripts/fit-model.R -i test-ising-123456.5.3.l1.counts -m genmatrices/genmatrix-5-complete.RData -j 54321
# 
# Rscript ../scripts/count-seq.R -i test-ising-123456.RData -w 6 -s 2 -l 2
# Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 6
# Rscript ../scripts/fit-model.R -i test-ising-123456.6.2.l2.counts -m genmatrices/genmatrix-6-complete.RData -j 54321
 
# ###
# # and for the constrained model
Rscript ../scripts/sim-seq.R -c ising-model-constrained.json -t .1 -s 100000 -o simseqs/constrained-ising-123.RData

# generator matrices
Rscript ../scripts/make-genmat.R -c ising-model-constrained.json -w 4 -o genmatrices/genmatrix-4-ising-model-constrained.RData
Rscript ../scripts/make-genmat.R -c ising-model-constrained.json -w 7 -o genmatrices/genmatrix-7-ising-model-constrained.RData
../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd simseqs/constrained-ising-123.RData genmatrices/genmatrix-7-ising-model-constrained.RData

# with shortish windows
Rscript ../scripts/count-seq.R -i simseqs/constrained-ising-123.RData -w 4 -s 2 -l 1 -o simseqs/constrained-ising-123.4.2.l1.counts
Rscript ../scripts/fit-model.R -i simseqs/constrained-ising-123.4.2.l1.counts -m genmatrices/genmatrix-4-ising-model-constrained.RData -o simseqs/constrained-ising-123-genmatrix-4-ising-model-constrained.RData
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/constrained-ising-123-genmatrix-4-ising-model-constrained.RData simseqs/constrained-ising-123.RData

# with longer windows
Rscript ../scripts/count-seq.R -i simseqs/constrained-ising-123.RData -w 7 -s 3 -l 2 -o simseqs/constrained-ising-123.7.3.l2.counts
Rscript ../scripts/fit-model.R -i simseqs/constrained-ising-123.7.3.l2.counts -m genmatrices/genmatrix-7-ising-model-constrained.RData -o simseqs/constrained-ising-123-genmatrix-7-ising-model-constrained.RData
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/constrained-ising-123-genmatrix-7-ising-model-constrained.RData simseqs/constrained-ising-123.RData

# Rscript ../scripts/mcmc-model.R -i constrained-ising-123-genmatrix-9-ising-model-constrained.RData -c ising-model-constrained.json -b 3 -o constrained-ising-123-genmatrix-9-mcmc-1.RData
# Rscript ../scripts/mcmc-model.R -i constrained-ising-123-genmatrix-9-ising-model-constrained.RData -c ising-model-constrained.json -b 100 -o constrained-ising-123-genmatrix-9-mcmc-2.RData
