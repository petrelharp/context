#!/bin/bash
# USAGE:
#    submit-job.sh (r script) (options)
# e.g.
#    submit-job.sh bcells-inference.R -u 02-C-M_out_of_frame -w 2 -l 2 -r 2 -k 2


read -d '' PBS_HEAD <<'EOF'
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=120:00:00
#PBS -l pmem=4gb
#PBS -l mem=500mb
#PBS -l vmem=5gb
EOF

NAME=$(echo "$*"| sed -e 's/[ \.]//g')

read -d '' SETUP <<'EOF'
source /usr/usc/R/3.0.2/setup.sh 
EOF


( echo "$PBS_HEAD"; 
        echo "#PBS -N $NAME"; 
        echo "$SETUP"; 
        echo "echo 'script:'";
        echo "echo '$*'";
        echo "";
        echo "cd $PWD";
        echo "Rscript --vanilla $* --jobid \$PBS_JOBID;" 
)  | qsub -
