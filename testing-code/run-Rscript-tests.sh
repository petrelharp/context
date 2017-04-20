#!/bin/bash

YESGOOD=0

LONGERSCRIPTs="check-ising-short-sim.R check-ising-long-sim.R check-ising-inference.R"  # like a minute each

for SCRIPT in $LONGERSCRIPTS
do
    echo "-------------- Running $SCRIPT"
    results=$(Rscript $SCRIPT | grep "Error:")
    if [ "$results" ]
    then
        echo "$SCRIPT failed with:"
        echo $results
        YESGOOD=1
    fi
done

exit $YESGOOD
