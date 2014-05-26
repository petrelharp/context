#!/bin/bash
# USAGE:
#    submit-job.sh (r script) (options)
# e.g.
#    submit-job.sh bcells-inference.R -u 02-C-M_out_of_frame -w 2 -l 2 -r 2 -k 2

if [[ ! $# > 0 ]]
then
    echo "Example: submit-job.sh bcells-inference.R -u 02-C-M_out_of_frame -w 2 -l 2 -r 2 -k 2"
    exit
fi

read -d '' PBS_HEAD <<'EOF'
#PBS -S /bin/bash
#PBS -q cmb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=120:00:00
#PBS -l pmem=4gb
#PBS -l mem=500mb
#PBS -l vmem=5gb
#PBS -j oe
EOF

NAME=$(echo "$*"| sed -e 's/[ \.]//g')


( echo "$PBS_HEAD"; 
        echo "#PBS -N $NAME"; 
        echo "#PBS -o qsub-logs/$NAME.o\$PBS_JOBID"; 
        echo "";
        echo "cd $PWD";
        echo "source /usr/usc/R/3.0.2/setup.sh";
        echo "echo 'script:'";
        echo "echo '$*'";
        echo 'git log -1 --format="%H"';
        echo "";
        echo "Rscript $* --jobid \$PBS_JOBID" 
)  | qsub -
