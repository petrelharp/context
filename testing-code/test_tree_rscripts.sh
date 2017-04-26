#!/bin/bash
set -eu
set -o pipefail

Rscript ../scripts/sim-seq.R -c fwds-rev-model.json -s 100 -o TEMP-tree-test.RData 2>/dev/null
Rscript -e 'load("TEMP-tree-test.RData");suppressMessages({require(Biostrings)})' \
    -e 'x <- lapply(simseqs,function (z) strsplit(as.character(z$finalseq),"")[[1]])' \
    -e 'stopifnot(unique(x$sp1)=="O")' \
    -e 'stopifnot(unique(x$sp2)=="O")' \
    -e 'stopifnot(unique(x$sp3)=="O")' \
    -e 'stopifnot(unique(x$sp4)=="X")' \
    -e 'stopifnot(unique(x$an5)=="X")' \
    -e 'stopifnot(all(x$an6==x$root))'

Rscript ../scripts/sim-seq.R -c fast-slow-tree.json -s 100 -o TEMP-tree-test.RData 2>/dev/null
Rscript -e 'load("TEMP-tree-test.RData");suppressMessages({require(Biostrings)})' \
    -e 'x <- lapply(simseqs,function (z) strsplit(as.character(z$finalseq),"")[[1]])' \
    -e 'stopifnot(all(x$sp1==x$an5))' \
    -e 'stopifnot(all(x$sp3==x$an6))' \
    -e 'stopifnot(all(x$an5==x$root))' 

rm TEMP-tree-test.RData
rm TEMP-tree-test.Rout
