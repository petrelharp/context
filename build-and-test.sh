set -eu

cd nestly && scons simple && scons -j 6 seed
cd ../cpg && ./minimal-workflow.sh
cd ../tasep && ./workflow.sh
cd ../ising && ./minimal-workflow.sh
