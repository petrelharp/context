#!/bin/bash

USAGE="Split a file up into chunks of N columns.
Usage:
   $0 N file
"

if [ $# -lt 2 ]
then
    echo "$USAGE"
    exit 0
fi

N=$1

shift

FILE=$1

shift

NCOL=$(head -n 1 $FILE | wc -w)

START=1
while [ $START -le $NCOL ]
do
    END=$((START+$N))
    cut $* -f $START-$END $FILE | column -t
    START=$((START+$N+1))
    echo ""
    echo "..........................................................."
    echo ""
done


