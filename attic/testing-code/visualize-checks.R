#!/usr/bin/R

library(contextual)
library(contextutils)
library(simcontext)

# Test suite?
require(parallel)
numcores <- detectCores()


###
# Simulate a bunch of short seqs
nsamples <- 10000
nuniqs <- 10
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
inseqs <- rep( sample(getpatterns(seqlen,bases),nuniqs,replace=TRUE), each=nsamples/nuniqs )
many.seqs <- mclapply(inseqs, function (initseq) { simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, bases=bases, initseq=initseq, fixfn=null.fixfn ) }, mc.cores=numcores )

# check transition matrix
leftwin <- 0
rightwin <- 0
shortwin <- seqlen
longwin <- leftwin+shortwin+rightwin
subtransmatrix <- gettransmatrix(mutpats, mutrates*tlen, selpats, selcoef, Ne, shortwin, leftwin, rightwin)
allseqs <- data.frame( 
        initseq=factor(sapply(many.seqs,"[[","initseq"),levels=rownames(subtransmatrix)), 
        finalseq=factor(sapply(many.seqs,"[[","finalseq"),levels=colnames(subtransmatrix)) 
    )
all.counts <- table(allseqs)
init.counts <- rowSums(all.counts)
nonz <- ( init.counts>0 )
expec.counts <- (init.counts*subtransmatrix)
if (interactive()) {
    layout(t(1:2))
    plot( as.vector(expec.counts[nonz]), as.vector(all.counts[nonz]) )
    abline(0,1)
    plot( as.vector(expec.counts[nonz]), as.vector(abs(all.counts[nonz]-expec.counts[nonz])) )
    abline(0,-qnorm(1/nsamples))
}

bignums <- as.vector( expec.counts > 1 )
stopifnot( all( abs(all.counts[bignums]-expec.counts[bignums])/expec.counts[bignums] < -qnorm(.01/nsamples)) )

if (interactive()) {
    plot( abs(all.counts[bignums]-expec.counts[bignums])/expec.counts[bignums] )
}

#####
# Really simple case

seqlen <- 1e5
Ne <- 1e4
tlen <- 6e6
patlen <- 2
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  lapply( combn( bases, 2, simplify=FALSE ), rev),  # single-base rates
        # list( c("CG","TG"), c("CG","CA") ), # CpG rates
    NULL )
# mutrates <- rep(1e-8,length(mutpats))
mutrates <- runif( length(mutpats) )*1e-8
selpats <- c( "[CG]" )
selcoef <- rep( 1e-4, length(selpats) )
simseqs <- simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef )

# full instantaneous mutation, and transition matrix
#   RECALL THIS OMITS THE DIAGONAL
genmatrix <- makegenmatrix( mutpats, selpats, patlen=patlen )
genmatrix@x <- update(genmatrix,mutrates/(patlen-sapply(sapply(mutpats,"[",1),nchar)+1),selcoef,Ne) # see source of simseq for why rescale
# check: number of (infinitesimal) transitions matches expected?
# normalized prob of transitions
obs.table <- table(simseqs$ntrans[c("i","j")])
obs.p <- with( simseqs, sweep(as.matrix(obs.table),1,rowSums(obs.table),"/") )
exp.p <- sweep(genmatrix,1,max(rowSums(genmatrix)),"/")
diag(exp.p) <- 1-rowSums(exp.p)
pvals <- rep( as.vector(exp.p), as.vector(obs.table) )
stopifnot(all(pvals>0))
stopifnot(all(pvals*length(pvals)>1/100))
if (interactive()) {
    plot(as.vector(exp.p),as.vector(obs.p)-as.vector(exp.p), xlab="expected infinitesimal transition counts", ylab="residuals")
    abline(0,5); abline(0,-5)
}

if (interactive()) {
# size of window on either side of the focal site
leftwin <- 2
rightwin <- 2
shortwin <- 1
longwin <- leftwin+shortwin+rightwin
subtransmatrix <- gettransmatrix(mutpats, mutrates, selpats, selcoef, Ne, tlen, shortwin, leftwin, rightwin)
counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs=simseqs, leftwin=leftwin )
# averaged
expected <- (seqlen-longwin+1) * (1/nbases)^longwin * subtransmatrix
# accounting for initial sequence
in.expected <- rowSums(counts) * subtransmatrix
# indicator of counts where patterns have changed
changed <- whichchanged(subtransmatrix,leftwin=leftwin,shortwin=shortwin)

tmp <- arrayInd( which(!changed & abs(counts-in.expected)<.5), .dim=dim(subtransmatrix) )
tmp <- data.frame( iseq=factor(tmp[,1],levels=1:nrow(subtransmatrix)), fseq=factor(tmp[,2],levels=1:ncol(subtransmatrix)) )
levels(tmp$iseq) <- rownames(subtransmatrix)
levels(tmp$fseq) <- colnames(subtransmatrix)
table(tmp$fseq)

# should be mixture of poissons -- compare means in bins
usethese <- ( changed & as.vector(in.expected)>.1 )
exp.bins <- cut( as.vector(in.expected)[usethese], 40 )
exp.bin.vals <- tapply( as.vector(in.expected)[usethese], exp.bins, mean )

    # observed = expected?
    layout(matrix(1:4,2))
    plot( as.vector(counts), as.vector(expected), col=1+changed ); abline(0,1)
    plot( as.vector(counts), as.vector(in.expected), col=1+changed ); abline(0,1)
    plot( as.vector(counts)[changed], as.vector(expected)[changed], col=2 ); abline(0,1)
    # should be mixture of poissons -- compare means in bins
    usethese <- ( changed & as.vector(expected)>.1 )
    exp.bins <- cut( as.vector(expected)[usethese], 40 )
    exp.bin.vals <- tapply( as.vector(expected)[usethese], exp.bins, mean )
    points( tapply( as.vector(counts)[usethese], exp.bins, mean ), exp.bin.vals, pch=20 )
    plot( as.vector(counts)[changed], as.vector(in.expected)[changed], col=2 ); abline(0,1)
    # should be mixture of poissons -- compare means in bins
    usethese <- ( changed & as.vector(in.expected)>.1 )
    exp.bins <- cut( as.vector(in.expected)[usethese], 40 )
    exp.bin.vals <- tapply( as.vector(in.expected)[usethese], exp.bins, mean )
    points( tapply( as.vector(counts)[usethese], exp.bins, mean ), exp.bin.vals, pch=20 )
}


###
# visualize changes

if (interactive()) {
    nsteps <- 50
    seqlen <- 200
    Ne <- 1e4
    tlen <- 6e6/nsteps
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
    inseq <- simseq( seqlen, tlen, genmatrix, bases )
    many.seqs <- character(nsteps)
    many.seqs[1] <- inseq$initseq
    many.seqs[2] <- show.simseq(inseq)
    for (k in 3:nsteps) { 
        inseq <- simseq( seqlen, tlen, genmatrix, bases, initseq=inseq$finalseq )
        many.seqs[k] <- show.simseq( inseq )
    }
    many.seqs
}
