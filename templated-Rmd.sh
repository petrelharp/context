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
outfile=$(pwd)/$(dirname $RData)/$(basename $RData .RData)-$(basename $template .Rmd).html

R --vanilla --slave << EOF
library("rmarkdown")

source("$(dirname $0)/context-inference-fns.R")
# Want to be able to load e.g. generator matrices and results from MCMC runs.
#  But: what's stored in these files might overwrite each other; 
#  specifically, results of fit-model and mcmc-model are all called 'model'.
#  Here, we rename the MCMC runs 'mcmcX', where X is an integer.
for (rdata in ${RDataList}) { 
    mcmcnum <- 1
    env <- new.env()
    load(rdata,env=env) 
    in.env <- ls(env=env)
    if ("model"%in%in.env && with(env,class(model)=="contextMCMC")) {
        assign(paste("mcmc",mcmcnum,sep=''),get("model",env))
        mcmcnum <- mcmcnum+1
        in.env <- setdiff(in.env,"model")
    }
    for (x in in.env) { assign(x,get(x,env)) }
}
render("$template", output_file="$outfile")
EOF
