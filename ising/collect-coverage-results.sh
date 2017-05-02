#!/bin/bash

if [ $# -gt 0 ]
then
    NJOBS=$1
else
    NJOBS=1
fi

SIMDIRS=$(find simseqs/coverage -maxdepth 1 -mindepth 1 -type d)

for DIR in $SIMDIRS
do
    SIMFITS=$(find $DIR -name "ising-fit*.RData")
    for SIMFIT in $SIMFITS
    do
        OUTFILE=$(echo $SIMFIT | sed -e 's/RData$/json/')
        echo "Reading from ${SIMFIT}, writing to ${OUTFILE}."
        Rscript ../scripts/gather-results.R --fit $SIMFIT --sim $DIR/ising.RData --outfile $OUTFILE --json &
        while (( $(jobs 2>&1 | grep -c Running) >= $NJOBS )); do sleep 0.1; done
    done
done

