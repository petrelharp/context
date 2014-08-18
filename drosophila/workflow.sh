#!/bin/bash

# download pairwise alignments
get-data.sh

# count tuples
../tuple-counting/count-tuples.sh data counts 3 1 1
../tuple-counting/count-tuples.sh data counts 6 2 2

# make generator matrices
