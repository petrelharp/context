#!/bin/bash

SIMDIRS=$(find simseqs -type d -name "sim-*")

for DIR in $SIMDIRS
do
    SIMFITS=$(find $DIR -name "test-ising-fit*.RData")
    for SIMFIT in $SIMFITS
    do
        OUTFILE=$(echo $SIMFIT | sed -e 's/RData$/json/')
        echo "Reading from ${SIMFIT}, writing to ${OUTFILE}."
        Rscript ../scripts/gather-results.R --fit $SIMFIT --sim $DIR/test-ising.RData --outfile $OUTFILE --json 2>/dev/null &
        while (( $(jobs 2>&1 | grep -c Running) >= 8 )); do sleep 0.1; done
    done
done

