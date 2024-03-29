#!/bin/bash

set -eu
set -o pipefail

if [ $# -ne 0 ]
then
    BETA=$1
else
    echo "Usage:   $0 (value of beta)"
    exit 0
fi

MODELFILE=ising-model-just-beta.json
BASEDIR=beta-runs/beta_$BETA
mkdir -p $BASEDIR

MODELNAME=${MODELFILE%.json}

# precompute generator matrix:
#   the `complete.json` file is the same as ising-model.json but without model fitting stuff
for LONGWIN in 6
do
    GENMATFILE=genmatrices/${MODELNAME}_genmatrix-${LONGWIN}.RData
    if ! [ -f $GENMATFILE ]
    then
        Rscript ../scripts/make-genmat.R -c $MODELFILE -w $LONGWIN -o $GENMATFILE
    fi
done

# this takes ~700 secs to simulate
SEQLEN=1000000
TLEN=1.0

MODEL=$BASEDIR/$MODELFILE
echo "Copying configuration, updated to $MODEL "
cat $MODELFILE | sed -e "s/\"selcoef\" : \[[ 0-9.]*\]/\"selcoef\" : [ $BETA ]/" > $MODEL

# number of MCMC batches of length 100 each
MCMCITER=10000

NRUNS=1

for N in $(seq $NRUNS)
do
    (
        DIR=${BASEDIR}/$(printf "%05g" $RANDOM);
        echo "Simulation $N, in $DIR"
        mkdir -p $DIR
        # simulate up some sequence for testing;
        Rscript ../scripts/sim-seq.R -c $MODEL -t $TLEN -s $SEQLEN -d $DIR -o ising.RData;
        # and count the Tmers;
        Rscript ../scripts/count-seq.R -i $DIR/ising.RData -w 6 -s 2 -l 2;
        LONGWIN=6

        # fit the model;
        GENMATFILE=genmatrices/${MODELNAME}_genmatrix-${LONGWIN}.RData
        Rscript ../scripts/fit-model.R -c $MODEL -i $DIR/ising-6-root-2-tip-l2-shift0.counts -t $TLEN -m $GENMATFILE -o $DIR/ising-fit-6-2-2.RData;

        # and mcmc
        MCMCID=$RANDOM
        Rscript ../scripts/mcmc-model.R -i $DIR/ising-fit-6-2-2.RData -c $MODEL -b $MCMCITER -l 1 -j $MCMCID

        # gather results into json
        for RDATA in $DIR/ising-fit*.RData
        do
            Rscript ../scripts/gather-results.R --json -f $RDATA -s $DIR/ising.RData > ${RDATA%RData}json
        done

    ) ##  NOT in parallel &
done

wait;

Rscript ../scripts/collect-params-results.R $(ls -t $BASEDIR/*/*/ising-fit*.json) > beta-experiment_results.tsv

PLOTSCRIPT="
all.res <- read.table('beta-experiment_results.tsv', header=TRUE, check.names=FALSE)
all.res <- all.res[grepl('mcmc', all.res[['file']]),]

paramnames <- c('O->X|X->O', 'OX|XO')
paramlabels <- c(expression(lambda), expression(beta))
true.cols <- match(paste0('sim:', paramnames), colnames(all.res))
est.cols <- match(paste0('q50%.', paramnames), colnames(all.res))
q.cols <- outer(c('q2.5%', 'q25%', 'q75%', 'q97.5%'), paramnames, 
                function(x,y) { match(paste(x,y,sep='.'), colnames(all.res)) } )

coverages <- tapply( 1:nrow(all.res), all.res[['longwin']], function (kk) { res <- all.res[kk,]; 
    sapply( paramnames, function (pn) {
        k <- match(pn, paramnames)
        ut <- TRUE
        c('95%'=1 - mean((res[ut,q.cols[1,k]] > res[ut,true.cols[k]]) | (res[ut,q.cols[4,k]] < res[ut,true.cols[k]])),
          '50%'=1 - mean((res[ut,q.cols[2,k]] > res[ut,true.cols[k]]) | (res[ut,q.cols[3,k]] < res[ut,true.cols[k]])))
      } ) } )

xtable::xtable(as.data.frame(lapply(coverages,function(x)x['95%',])))

    beta <- res[,'sim:OX|XO']
    bvals <- sort(unique(beta))
    res <- all.res[match(bvals,beta),]

    pdf(file='beta-results.pdf', width=8, height=4, pointsize=10)
layout(t(1:2))
matplot(bvals, res[,c('q2.5%.OX|XO','q25%.OX|XO','q50%.OX|XO','q75%.OX|XO','q97.5%.OX|XO')]-bvals, xlab=expression(beta), ylab='error', type='l')
abline(v=1)
matplot(bvals, res[,c('q2.5%.OX|XO','q25%.OX|XO','q50%.OX|XO','q75%.OX|XO','q97.5%.OX|XO')]-bvals, ylim=c(-.1,.1), xlim=c(0,2), xlab=expression(beta), ylab='error', type='l')
abline(v=1)
dev.off()

#     pdf(file=sprintf('../writeups/writeup-plots/beta_results.pdf'),
#         width=5.5, height=3, pointsize=10)
    xoffset <- 0
    xshift <- nrow(res)/5
    par(mar=c(2,4,1,1)+.1)
    plot(0, type='n', 
          xlim=c(0, length(paramnames)*nrow(res) + (length(paramnames)-1)*xshift),
          ylim=c(0.85,1.15), xaxt='n',
          xlab='', ylab='relative parameter value')
    abline(h=1, col='red', lty=1)
    for (k in seq_along(paramnames)) {
        if (k>1) abline(v=xoffset - xshift/2, lty=2)
        xord <- order(res[['sim:OX|XO']],res[,est.cols[k]])
        axis(1, at=xoffset+(nrow(res)+xshift)/2, 
             labels=paramlabels[k], 
             mgp=c(3,0.2,0), tick=FALSE)
        segments(x0=xoffset+(1:nrow(res)), 
                 y0=res[xord,q.cols[1,k]]/res[xord,true.cols[k]], 
                 y1=res[xord,q.cols[4,k]]/res[xord,true.cols[k]],
                 col='grey', lwd=0.5)
        segments(x0=xoffset+(1:nrow(res)), 
                 y0=res[xord,q.cols[2,k]]/res[xord,true.cols[k]], 
                 y1=res[xord,q.cols[3,k]]/res[xord,true.cols[k]],
                 col='black', lwd=1)
    #     points(x=xoffset+(1:nrow(res)), y=res[xord,est.cols[k]]/res[xord,true.cols[k]], 
    #             cex=0.5, pch=20)
        xoffset <- xoffset + nrow(res) + xshift
        abline(v=xoffset - xshift/2, lty=2)
    }
    abline(h=1, col='red', lty=3)
#     dev.off()

"
# echo "$PLOTSCRIPT" | littler
