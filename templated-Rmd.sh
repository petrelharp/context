#!/bin/sh

set -e
set -o

test "$#" -eq 2 || {
  echo "Usage: $0 TEMPLATE RDATA" >&2
  exit 1
}

template=$1
RData=$2
output=$(basename $RData .RData).html

R --vanilla --slave << EOF
library("rmarkdown")

source("$(dirname $0)/context-inference-fns.R")
load("$RData")
render("$template",output_file="$output")
EOF

