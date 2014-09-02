#!/bin/bash

set -eu 
set -o pipefail

Rscript ../sim-seq.R -c tree-cpg-model.json -s 10000 -o sim-tree-cpg-model-111.RData
