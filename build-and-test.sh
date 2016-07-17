set -eu

cd nestly && scons simple && scons -j 6 seed
cd ../json-cpg && ./minimal-workflow.sh
cd ../json-tasep && ./workflow.sh
cd ../json-ising && ./minimal-workflow.sh
cd ../testing-code && ./run-Rscript-tests.sh
