#!/bin/bash

DIRS="RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx"
TYPES="CTCF_binding_site  enhancer  open_chromatin_region  promoter_flanking_region"

## Base model:
# get the correct base frequencies in there:
for d in $DIRS; do for t in $TYPES; do 
    CG=$(grep "$(basename $d)/$t" RegulatoryFeature-regions-from-axt/initfreqs.txt | cut -f 2 -d ' ')
    # JSON needs leading 0.'s ...
    A=0$(bc -l <<< "scale=4; (1-$CG)/2.0")
    C=0$(bc -l <<< "scale=4; $CG/2.0")
    G=0$(bc -l <<< "scale=4; $CG/2.0")
    T=0$(bc -l <<< "scale=4; (1-$CG)/2.0")
    cat base-model.json | sed -e 's_genmatrices_../../genmatrices_' | sed -e "s/0.25, 0.25, 0.25, 0.25/$A, $C, $G, $T/" > $d/$t/base-model.json; 
done; done

# collapse down to 6-3-2
collapse () {
    DIR=$1
    COLLAPSE=" library(contextual); f <- file.path(dir, \"7-5-1/total.counts.gz\"); g <- file.path(dirname(dirname(f)), \"6-3-2\", \"total.counts.gz\"); c <- read.counts(f, leftwin=1); d <- projectcounts(c, new.leftwin=2, new.shortwin=3, new.longwin=6); dir.create(dirname(g), recursive=TRUE, showWarnings=FALSE); write.counts(d, file=g); "
    Rscript -e "dir=\"$DIR\"" -e "$COLLAPSE"
}

for d in $DIRS; do for t in $TYPES; do 
    collapse $d/$t
done; done

# create genmatrix
LONGWIN=6
../scripts/make-genmat.R -c base-model.json -w $LONGWIN -o genmatrices/base-model-genmatrix-${LONGWIN}.RData

# fit model
TMER="6-3-2"
for d in $DIRS; do for t in $TYPES; do 
    ../scripts/fit-tree-model.R -i $d/$t/$TMER/total.counts.gz -c $d/$t/base-model.json -o $d/$t/$TMER/base-model-fit.RData 
done; done

