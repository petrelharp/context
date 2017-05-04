---
title: Methods for getting axt alignments of regulatory regions
author: Jessica Crisci
date: July 30, 2015
---

How files in this directory were obtained
====

* Human-chimp axt alignments (hg38-panTro4) were downloaded from [UCSC] (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro4/).
* Regulatory features in hg38 for multiple cell types were downloaded from the [Ensembl Regulatory Build] (ftp://ftp.ensembl.org/pub/release-81/regulation/homo_sapiens/).
* The knownGenes table was downloaded from the [hg38 Annotations Database] (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/).

Filtering regulatory regions for length
----

Next, I made histograms of the region lengths for each of the six regualtory feature categories in the Ensembl Regulatory Build file

![histograms] (hist_and_filter-regulatory-region-length.jpg)

* I chose a region length filter of < 1000bp for CTCF_binding_site, TFBS, enhancer, and open_chromatin_region.
* For promoter and promoter_flanking_regions I chose a filter of < 3000bp.

Sorting regulatory regions by overlap with knownGenes
----

After filtering for length, I sorted the regulary regions in each category by those that overlapped with transcription start and end regions from the UCSC
knownGene table and those that shared no overlap (see subdirectories `overlap-knownGeneTx` and `noOverlap-knownGeneTx`). Within each of these subdirectories
is a series of files, per chromosome, of axt alignements that overlap with each regulatory feature, filtered for length, and sorted by overlap with transcribed
regions of the genome. 

Data notes
------

- Lower case bases are (soft) repeat-masked locations.

Counting tuples
------

In the files, `reference` is *human* and `derived` is *chimp*.

*3 August 2015, plr*: 
Syntax:
```
Usage:
    count-tuples.sh (indir) (outdir) (longwin) (leftwin) (rightwin)
```
Note that `count-tuples.sh` calls `count-paired-tuples.py` with `--strict`, which ignores anything that isn't A, C, G, or T (case-sensitive),
so this is ignoring the repeat-masked bases.
I have genmatrices for lengths 5, 7, 8, and 9, so,
to count paired tuples in each .axt file:
```{#sh}
COUNTSCRIPT="/home/peter/projects/context/tuple-counting/count-tuples.sh"
for dir in RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx;
do
    # leftwin+shortwin+rightwin
    # 0+5+0 = 5
    ( LW=5; SW=5; LT=0; RT=$(($LW - $SW - $LT)); OUTDIR=${dir}-${LW}-${SW}-${LT}; mkdir -p $OUTDIR; if [ ! -e $OUTDIR/count-tuples.log ]; then ${COUNTSCRIPT} $dir $OUTDIR $LW $LT $RT &> $OUTDIR/count-tuples.log; fi; ) &
    # 1+5+1 = 7
    ( LW=7; SW=5; LT=1; RT=$(($LW - $SW - $LT)); OUTDIR=${dir}-${LW}-${SW}-${LT}; mkdir -p $OUTDIR; if [ ! -e $OUTDIR/count-tuples.log ]; then ${COUNTSCRIPT} $dir $OUTDIR $LW $LT $RT &> $OUTDIR/count-tuples.log; fi; ) &
    # 1+6+1 = 8
    ( LW=8; SW=6; LT=1; RT=$(($LW - $SW - $LT)); OUTDIR=${dir}-${LW}-${SW}-${LT}; mkdir -p $OUTDIR; if [ ! -e $OUTDIR/count-tuples.log ]; then ${COUNTSCRIPT} $dir $OUTDIR $LW $LT $RT &> $OUTDIR/count-tuples.log; fi; ) &
    # 2+5+2 = 9
    ( LW=9; SW=5; LT=2; RT=$(($LW - $SW - $LT)); OUTDIR=${dir}-${LW}-${SW}-${LT}; mkdir -p $OUTDIR; if [ ! -e $OUTDIR/count-tuples.log ]; then ${COUNTSCRIPT} $dir $OUTDIR $LW $LT $RT &> $OUTDIR/count-tuples.log; fi; ) &
done
```
These finished in about 10 minutes.

Hm, these are split up by chromosome, and we actually want to agglomerate these across chromosomes.
To do this, use `sum-counts.py`:
```{#sh}
SUMSCRIPT="/home/peter/projects/context/tuple-counting/sum-counts.py"
for dir in RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx
do
    for tuples in $(ls -d ${dir}-* | sed -e 's/.*Tx-//g')
    do
        ( for type in CTCF_binding_site enhancer open_chromatin_region promoter_flanking_region
        do
            mkdir -p $dir/$type/$tuples
            $SUMSCRIPT -i ${dir}-${tuples}/${type}.*.counts.gz -o $dir/$type/$tuples/total.counts.gz
        done ) &
    done
done
```
This takes less than a minute.
These will then be moved to the `context` project:
```{#sh}
CDIR=/home/peter/projects/context/shape
mkdir -p $CDIR/RegulatoryFeature-regions-from-axt
cp -r $(find RegulatoryFeature-regions-from-axt -mindepth 2 -maxdepth 2 -type 'd') $CDIR/RegulatoryFeature-regions-from-axt
```
