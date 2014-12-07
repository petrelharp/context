#!/bin/bash

YESGOOD=0

for SCRIPT in sim-short-ising.R check-projection.R tree-tests.R check-ising.R
do
    results=$(Rscript $SCRIPT | grep "Error:")
    if [ "$results" ]
    then
        echo "$SCRIPT failed with:"
        echo $results
        YESGOOD=1
    fi
done

exit $YESGOOD
