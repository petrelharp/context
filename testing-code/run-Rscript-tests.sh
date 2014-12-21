#!/bin/bash

YESGOOD=0

QUICKSCRIPTS="tree-tests.R tree-counts-test.R check-projection.R sim-short-ising.R check-ising-genmats.R"
LONGERSCRIPTs="check-ising-short-sim.R check-ising-long-sim.R check-ising-inference.R"  # like a minute each

for SCRIPT in $QUICKSCRIPTS $LONGERSCRIPTS
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
