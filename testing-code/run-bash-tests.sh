#!/bin/bash

set -eu

./unit-tests.sh
./unit-test-treeness.sh
./test-sim-and-save.sh


