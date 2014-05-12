#!/bin/bash

basename=$1
fasta="${basename}.fasta"

if [ ! -d "$basename" ]
then
    mkdir $basename
fi

w=5; lr=2; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=5; lr=1; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=6; lr=1; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=6; lr=2; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=7; lr=3; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=7; lr=2; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=8; lr=3; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=8; lr=2; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
w=9; lr=3; python ../tuple-counting/count-paired-tuples.py -w ${w} -l ${lr} -r ${lr} -i $fasta > $basename/tuples.${w}.${lr}.counts
