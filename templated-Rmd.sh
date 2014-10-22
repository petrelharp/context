#!/bin/sh

set -e
set -u

usage="Usage: $0 TEMPLATE RDATA [RDATA2]"

test "$#" -ge 2 || {
  echo $usage >&2
  exit 2
}

template=$1
RData=$2
shift
RDataList="c(\""$(echo $@| sed -e 's/ /", "/g')"\")"
outfile=$(pwd)/$(dirname $RData)/$(basename $RData .RData).html

R --vanilla --slave << EOF
library("rmarkdown")

source("$(dirname $0)/context-inference-fns.R")
for (rdata in ${RDataList}) { load(rdata) }
render("$template", output_file="$outfile")
EOF
