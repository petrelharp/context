#!/bin/bash
# Run is-coding.py on everything.

SCRIPTDIR=/home/peter/projects/codons

codons="GCT GCC GCA GCG CGT CGC CGA CGG AGA AGG AAT AAC GAT GAC TGT TGC CAA CAG GAA GAG GGT GGC GGA GGG CAT CAC ATT ATC ATA TTA TTG CTT CTC CTA CTG AAA AAG ATG TTT TTC CCT CCC CCA CCG TCT TCC TCA TCG AGT AGC ACT ACC ACA ACG TGG TAT TAC GTT GTC GTA GTG"

for chrom in 2L 2R 3L 3R 4 X
do
    for codon in $codons
    do
        python $SCRIPTDIR/is-coding.py -i ${chrom}.raw.${codon}.gz -c CDS-dmel-${chrom}-r5.50.CDS.starts.ends.gz -o ${chrom}.raw.${codon}
        # check:
        nnc=$( zcat ${chrom}.raw.${codon}.noncoding.gz | wc -l )
        noc=$( zcat ${chrom}.raw.${codon}.nonorfcoding.gz | wc -l )
        oc=$( zcat ${chrom}.raw.${codon}.orfcoding.gz | wc -l )
        echo $chrom $codon " : " $(zcat ${chrom}.raw.${codon}.gz | wc -l) "=" $nnc "+" $noc "+" $oc "=" $(($nnc+$noc+$oc))
    done
done

