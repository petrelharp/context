#!/bin/bash

usage="Usage: $0 CONFIG OUTFILE longwin=INT shortwin=INT leftwin=INT tlen=NUMERIC ncounts=INT [VAR=VALUE [VAR=VALUE]]"

test "$#" -ge 2 || {
  echo $usage >&2
  exit 2
}

pandoc -v >/dev/null || {
    echo "No working pandoc found.  Exiting."
    exit 3
}

template="$(dirname $0)/likelihood-surface.Rmd"
config="$(pwd)/$1"
outfile="$(pwd)/$2"
shift
shift

mkdir -p $(dirname $outfile)

# read does NOT play well with set -eu 
read -r -d '' SCRIPT << EOF
library("rmarkdown")
modelfile <- "${config}"
$(for x in "$@"; do echo $x; done)
render("$template", output_file="$outfile")
EOF

echo "$SCRIPT"

echo "$SCRIPT" | R --slave 
