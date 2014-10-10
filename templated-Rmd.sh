#!/bin/sh

set -e
set -u

usage="Usage: $0 [-o OUTFILE] TEMPLATE RDATA"

test "$#" -eq 2 || {
  echo $usage >&2
  exit 2
}

template=$1
RData=$2
outfile=$(pwd)/$(dirname $RData)/$(basename $RData .RData).html

R --vanilla --slave << EOF
library("rmarkdown")

source("$(dirname $0)/context-inference-fns.R")
load("$RData")
render("$template", output_file="$outfile")
EOF
