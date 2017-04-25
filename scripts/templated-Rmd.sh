#!/bin/sh

set -e
set -u

usage="Usage: $0 TEMPLATE RDATA [RDATA2]"

test "$#" -ge 2 || {
  echo $usage >&2
  exit 2
}

pandoc -v >/dev/null || {
    echo "No working pandoc found.  Exiting."
    exit 3
}

template=$1
RData=$2
shift
RDataList="c(\""$(echo $@| sed -e 's/ /", "/g')"\")"
outfile=$(pwd)/$(dirname $RData)/$(basename $RData .RData)-$(basename $template .Rmd).html

R --slave << EOF
library("rmarkdown")

library(contextual)
library(contextutils)
for (rdata in ${RDataList}) { cat("Loading ", rdata, "\n"); load(rdata) }
render("$template", output_file="$outfile")
EOF
