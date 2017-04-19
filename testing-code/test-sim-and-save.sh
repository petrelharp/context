#!/bin/bash

# test that we can simulate, saving the generator matrix, and the resimulate using it
rm -f minimal/tasep-genmat.RData

Rscript ../sim-seq.R -c ../json-tasep/tasep-model.json -t 0.1 -s 100 -d minimal -o sim-tasep-12345.RData -m minimal/tasep-genmat.RData && ls minimal/tasep-genmat.RData && ls minimal/sim-tasep-12345.RData
Rscript ../sim-seq.R -c ../json-tasep/tasep-model.json -t 0.1 -s 100 -d minimal -o sim-tasep-23456.RData -m minimal/tasep-genmat.RData && ls minimal/sim-tasep-23456.RData

