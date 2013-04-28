#!/usr/bin/R
source("codons.R")
source("codon-inference-fns.R")
source("sim-context-fns.R")

# Test suite?

###
# Simulate a bunch of short seqs
nsamples <- 10000
seqlen <- 4
Ne <- 1e4
tlen <- 6e7
patlen <- 2
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  lapply( combn( bases, 2, simplify=FALSE ), rev),  # single-base rates
        list( c("CG","TG"), c("CG","CA") ), # CpG rates
    NULL )
mutrates <- runif( length(mutpats) )*1e-8
selpats <- c( "[GC]", "[AT]" )
selcoef <- runif( length(selpats) )*1e-4
genmatrix <- makegenmatrix( mutpats, selpats, patlen=patlen )
genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)
inseqs <- rep( sample(getpatterns(seqlen),nsamples/100,replace=TRUE), each=100 )
many.seqs <- lapply(inseqs, function (initseq) { simseq( seqlen, tlen, genmatrix, bases, initseq=initseq ) }  )

# check transition matrix
lwin <- 0
rwin <- 0
win <- seqlen
winlen <- lwin+win+rwin
subtransmatrix <- gettransmatrix(mutpats, mutrates, selpats, selcoef, Ne, tlen, win, lwin, rwin)
allseqs <- data.frame( 
        initseq=factor(sapply(many.seqs,"[[","initseq"),levels=rownames(subtransmatrix)), 
        finalseq=factor(sapply(many.seqs,"[[","finalseq"),levels=colnames(subtransmatrix)) 
    )
all.counts <- table(allseqs)
init.counts <- rowSums(all.counts)
nonz <- ( init.counts>0 )
expec.counts <- (init.counts*subtransmatrix)
layout(t(1:2))
plot( as.vector(expec.counts[nonz]), as.vector(all.counts[nonz]) )
abline(0,1)
plot( as.vector(expec.counts[nonz]), as.vector(abs(all.counts[nonz]-expec.counts[nonz])) )
abline(0,-qnorm(1/nsamples))
bignums <- as.vector( expec.counts > 1 )
stopifnot( all( abs(all.counts[bignums]-expec.counts[bignums])/expec.counts[bignums] < -qnorm(.01/nsamples)) )

