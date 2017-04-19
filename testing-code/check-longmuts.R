#!/usr/bin/Rscript

# a contrived example to test longer mutpats

library(contextual)
library(contextutils)
library(simcontext)

bases <- c("O","X")

mutpats <- list( 
    list( c("O","X") ),
    list( c("OX","OO"), c("XO","OO") )
    ) 
mutrates <- c(3,5)
selpats <- list(
        c( "X" )
    )
selcoef <- c(1)

fixfn <- function (ds,...) { 1/(1+exp(-ds)) }


###### 
# Can make correct generator matrices?
patlen <- 3

pats <- c("OOO", "XOO", "OXO", "XXO", "OOX", "XOX", "OXX", "XXX")
nA <- rbind(   # number of O->X:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  1,  1,  0,  1,  0,  0,  0 ), # OOO
        c(  -1,  0,  0,  1,  0,  1,  0,  0 ), # XOO
        c(  -1,  0,  0,  1,  0,  0,  1,  0 ), # OXO
        c(   0, -1, -1,  0,  0,  0,  0,  1 ), # XXO
        c(  -1,  0,  0,  0,  0,  1,  1,  0 ), # OOX
        c(   0, -1,  0,  0, -1,  0,  0,  1 ), # XOX
        c(   0,  0, -1,  0, -1,  0,  0,  1 ), # OXX
        c(   0,  0,  0, -1,  0, -1, -1,  0 )  # XXX
    )
nB <- rbind(   # number of OX->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # XXO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
nC <- rbind(   # number of XO->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XXO
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # XOX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
dimnames(nA) <- dimnames(nB) <- dimnames(nC) <- list(pats,pats)

checkit <- function (x,y) { stopifnot( isTRUE(all.equal(unlist(dimnames(x)),unlist(dimnames(y)))) && isTRUE(all.equal(as.vector(x),as.vector(y))) ) }

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=null.fixfn ),
        3 * (nA>0) + 5 * ( nB + nC )
    )


## wrapped

nA.wrap <- nA
nB.wrap <- rbind(   # number of OX->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  0,  1,  0,  0,  0,  0,  0 ), # XXO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
nC.wrap <- rbind(   # number of XO->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XXO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # XOX
        c(   0,  0,  1,  0,  0,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
dimnames(nA.wrap) <- dimnames(nB.wrap) <- dimnames(nC.wrap) <- list(pats,pats)

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=ising.fixfn ),
        ( 3 * (nA.wrap>0) + 5 * ( nB.wrap + nC.wrap ) ) * ising.fixfn( 1 * nA )
    )

cat("\n\n  Generator matrices all good!\n")

# Short sequences:
require(parallel)
numcores <- getcores()

check.sim <- function (tlen, seqlen, mutrates, selcoef, initseq) {
    if (missing(initseq)) { initseq <- rinitseq(seqlen,bases) }
    simseqs <- mclapply( 1:1000, function (k) 
            simseq( seqlen, tlen, initseq=initseq, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, bases=bases, fixfn=fixfn, count.trans=TRUE ) 
        , mc.cores=numcores)
    # note that simseq simulates WRAPPED sequences
    full.genmatrix <- makegenmatrix( mutpats, selpats, mutrates=mutrates, selcoef=selcoef, patlen=seqlen, boundary="wrap", bases=bases, fixfn=fixfn )
    require(expm)
    full.transmatrix <- expm( tlen*(as.matrix(full.genmatrix)-diag(rowSums(full.genmatrix))) )
    stopifnot( all( abs( rowSums(full.transmatrix) - rep(1.0,nrow(full.transmatrix)) ) < 1e-8 ) )
    expected.freqs <- full.transmatrix[match(as.character(initseq),rownames(full.transmatrix)),]
    fseqs <- as.numeric(table( factor( sapply( simseqs, function (x) as.character(x$finalseq) ), levels=names(expected.freqs) ) ))
    # plot( length(simseqs)*expected.freqs, fseqs ); abline(0,1)
    resids <- data.frame( observed=fseqs, expected=expected.freqs*length(simseqs) )
    resids$z <- (resids$observed-resids$expected)/sqrt(resids$expected)
    resids$p <- pbinom( resids$observed, size=length(simseqs), prob=expected.freqs )
    return( list( initseq=initseq, simseqs=simseqs, resids=resids ) )
}


# with all parameters
full.sims <- check.sim(tlen=2, seqlen=5, mutrates=mutrates, selcoef=selcoef)

cat("Beginning from ", as.character(full.sims$initseq), ":\n")
print(full.sims$resids)

if ( with(full.sims, min( pbeta( min(resids$p), 1,nrow(resids) ) , pbeta( min(1-resids$p), 1,nrow(resids) ) ) < .01 ) ) {
    stop("some high residuals there.")
}
