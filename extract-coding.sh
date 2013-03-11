#!/bin/bash
# Extract the desired information from a GFF file

for chrom in 2L 2R 3L 3R 4 X
do
    zcat CDS-dmel-all-r5.50.gff.gz | grep $chrom | grep 'FlyBase	CDS' | cut -f 4-5 | gzip -c > CDS-dmel-${chrom}-r5.50.CDS.starts.ends.gz
done
