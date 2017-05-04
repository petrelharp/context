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

Getting average base frequencies
--------------------------------

TLDR: here are GC contents, which should suffice (C/G and A/T are nearly equal):
```
  overlap-knownGeneTx/enhancer      overlap-knownGeneTx/open_chromatin_region          overlap-knownGeneTx/CTCF_binding_site   overlap-knownGeneTx/promoter_flanking_region 
                     0.4251248                                      0.4578896                                      0.4478951                                      0.4393010 
noOverlap-knownGeneTx/enhancer    noOverlap-knownGeneTx/open_chromatin_region        noOverlap-knownGeneTx/CTCF_binding_site noOverlap-knownGeneTx/promoter_flanking_region 
                     0.4341862                                      0.4445058                                      0.4449049                                      0.4404174 
```

Doing
```

for DIR in overlap-knownGeneTx noOverlap-knownGeneTx
do
    for TYPE in enhancer open_chromatin_region CTCF_binding_site promoter_flanking_region
    do
        ( echo $DIR $TYPE; zcat $DIR/${TYPE}*.axt.gz | awk 'BEGIN {FS=""} {for (i=1; i<=NF; i++){a[$i]++}}END{for (i in a){print i, a[i]}}' | grep "[ACGT]" ) >> basefreqs.txt
    done
done

```
gets after parsing
```
                                         names        A        C        G        T        pA        pC        pG        pT
                  overlap-knownGeneTx enhancer 10406409  7692217  7715027 10428045 0.2871391 0.2122477 0.2128771 0.2877361
     overlap-knownGeneTx open_chromatin_region  5228168  4417041  4427844  5243577 0.2706563 0.2286652 0.2292244 0.2714540
         overlap-knownGeneTx CTCF_binding_site 10940737  8867546  8897490 10957615 0.2758397 0.2235701 0.2243250 0.2762652
  overlap-knownGeneTx promoter_flanking_region 20496007 16050051 16094793 20531859 0.2801045 0.2193448 0.2199562 0.2805945
                noOverlap-knownGeneTx enhancer  7190366  5510859  5512503  7174828 0.2832129 0.2170607 0.2171255 0.2826009
   noOverlap-knownGeneTx open_chromatin_region  4190742  3358080  3355372  4198988 0.2774741 0.2223425 0.2221632 0.2780201
       noOverlap-knownGeneTx CTCF_binding_site  8417303  6740981  6742683  8405881 0.2777360 0.2224243 0.2224805 0.2773591
noOverlap-knownGeneTx promoter_flanking_region 14310102 11266049 11262795 14314450 0.2797488 0.2202405 0.2201769 0.2798338
```
which is in R
```
structure(list(names = structure(c(6L, 7L, 5L, 8L, 2L, 3L, 1L, 
4L), .Label = c("noOverlap-knownGeneTx CTCF_binding_site", "noOverlap-knownGeneTx enhancer", 
"noOverlap-knownGeneTx open_chromatin_region", "noOverlap-knownGeneTx promoter_flanking_region", 
"overlap-knownGeneTx CTCF_binding_site", "overlap-knownGeneTx enhancer", 
"overlap-knownGeneTx open_chromatin_region", "overlap-knownGeneTx promoter_flanking_region"
), class = "factor"), A = c(10406409, 5228168, 10940737, 20496007, 
7190366, 4190742, 8417303, 14310102), C = c(7692217, 4417041, 
8867546, 16050051, 5510859, 3358080, 6740981, 11266049), G = c(7715027, 
4427844, 8897490, 16094793, 5512503, 3355372, 6742683, 11262795
), T = c(10428045, 5243577, 10957615, 20531859, 7174828, 4198988, 
8405881, 14314450), pA = c(0.287139112521715, 0.270656320486544, 
0.275839698817459, 0.280104522574058, 0.283212877487006, 0.277474111084671, 
0.277736008706679, 0.279748816676805), pC = c(0.212247697087906, 
0.228665193572084, 0.223570057563678, 0.219344765809009, 0.217060747311726, 
0.222342545981747, 0.222424349052226, 0.220240488791555), pG = c(0.212877081357832, 
0.229224449893346, 0.224325009451283, 0.219956222615549, 0.217125499038931, 
0.222163242749714, 0.222480506343783, 0.220176875254068), pT = c(0.287736098762701, 
0.271454016680615, 0.276265224747247, 0.2805944839017, 0.28260086147721, 
0.278020075456622, 0.277359123573951, 0.279833811981381)), .Names = c("names", 
"A", "C", "G", "T", "pA", "pC", "pG", "pT"), row.names = c(NA, 
-8L), class = "data.frame")
```
