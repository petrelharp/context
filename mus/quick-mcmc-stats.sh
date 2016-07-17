#!/bin/bash

for x in $(find . -name "win-?-?-?-mcmc-?.RData"); do (echo -n $x " " ; Rscript -e "load(\"$x\"); cat(c( lwin, win, rwin, mrun[['accept']], mrun[['nbatch']] * mrun[['blen']], mrun[['final']] )); cat('\n')" | column -t;) done
