#!/bin/bash

set -eu 
set -o pipefail

Rscript ../sim-seq.R -c tree-cpg-model.json -s 100 
