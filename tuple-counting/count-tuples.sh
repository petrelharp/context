#!/bin/bash

if [ $# -lt 5 ]
then
    echo "
Usage:
    count-tuples.sh (indir) (outdir) (longwin) (leftwin) (rightwin) 
e.g.
    count-tuples.sh data counts 4 1 1
    
This will then run count-paired-tuples.py on all .axt.gz files in the directory 'data', and place resulting files in the directory 'counts'.
"
    exit 1
fi

INDIR=$1
OUTDIR=$2
LONGWIN=$3
LEFTWIN=$4
RIGHTWIN=$5

for x in ${INDIR}/*.axt.gz
do
    echo "counting: " $x
    OUTNAME=$(echo $x | sed -e "s_${INDIR}/_${OUTDIR}/_" | sed -e "s/.axt.gz$//").${LONGWIN}.l${LEFTWIN}.r${RIGHTWIN}.counts.gz
    echo "   to: " $OUTNAME
    python ../tuple-counting/count-paired-tuples.py -i $x -f axt -o $OUTNAME -w $LONGWIN -l $LEFTWIN -r $RIGHTWIN --strict
done
