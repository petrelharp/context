library(testthat)

context("a contrived example to test longer mutpats")

library(contextual)
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


# Short sequences:

check.sim <- function (tlen, seqlen, mutrates, selcoef, initseq) {
    if (missing(initseq)) { initseq <- rinitseq(seqlen,bases) }
    simseqs <- lapply( 1:100, function (k) 
            simseq( seqlen, tlen, initseq=list(initseq), mutpats=mutpats, selpats=selpats, mutrates=mutrates, 
                   selcoef=selcoef, bases=bases, fixfn=fixfn, count.trans=TRUE, quiet=TRUE ) 
        )
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

# cat("Beginning from ", as.character(full.sims$initseq), ":\n")
# print(full.sims$resids)

context("This should fail once out of a thousand times.")

expect_true( with(full.sims, min( pbeta( min(resids$p), 1,nrow(resids) ) , pbeta( min(1-resids$p), 1,nrow(resids) ) ) > .001 ) )

