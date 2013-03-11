#!/bin/bash

codons="GCT GCC GCA GCG CGT CGC CGA CGG AGA AGG AAT AAC GAT GAC TGT TGC CAA CAG GAA GAG GGT GGC GGA GGG CAT CAC ATT ATC ATA TTA TTG CTT CTC CTA CTG AAA AAG ATG TTT TTC CCT CCC CCA CCG TCT TCC TCA TCG AGT AGC ACT ACC ACA ACG TGG TAT TAC GTT GTC GTA GTG"

for chrom in 2 3 4 X
do
    python get-positions.py -i ${chrom}.raw.gz -l $chrom.raw.patterns.log -o ${chrom}.raw --patterns $codons 
done
