#!/bin/sh

set -e
set -o

test "$#" -eq 2 || {
  echo "Usage: $0 TEMPLATE RDATA" >&2
  exit 1
}

template=$1
RData=$2

R --vanilla --slave << EOF
library("rmarkdown")

load("$RData")
render("$template")
EOF

