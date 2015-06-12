#!/usr/bin/Rscript

require(expm)
require(parallel)

# Run some tests!!
source("../sim-context-fns.R",chdir=TRUE)
source("../context-inference-fns.R",chdir=TRUE)
numcores<-getcores()

bases <- c("X","O")

patlen <- 2
mutpats <- list( 
    list( c("O","X") ),
    list( c("X","O") )
    ) 
selpats <- list(
        c("OX","XO"),
        c("X")
    )


fixfn <- function (ds,...) { 1/(1+exp(-ds)) }

# Short sequences:
tlen <- 2
seqlen <- 5
mutrates <- c(1,1)
selcoef <- c(-.5,.5)

initseq <- rinitseq(seqlen,bases)
simseqs <- mclapply( 1:1000, function (k) 
        simseq( seqlen, tlen, initseq=initseq, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, bases=bases, fixfn=fixfn, count.trans=TRUE ) 
    , mc.cores=numcores)
# note that simseq simulates WRAPPED sequences
full.genmatrix <- makegenmatrix( mutpats, selpats, murates=mutrates, selcoef=selcoef, patlen=seqlen, boundary="wrap", bases=bases, fixfn=fixfn )
require(expm)
full.transmatrix <- expm( tlen*(as.matrix(full.genmatrix)-diag(rowSums(full.genmatrix))) )
stopifnot( all( abs( rowSums(full.transmatrix) - rep(1.0,nrow(full.transmatrix)) ) < 1e-8 ) )
expected.freqs <- full.transmatrix[match(as.character(initseq),rownames(full.transmatrix)),]
fseqs <- as.numeric(table( factor( sapply( simseqs, function (x) as.character(x$finalseq) ), levels=names(expected.freqs) ) ))
# plot( length(simseqs)*expected.freqs, fseqs ); abline(0,1)
resids <- data.frame( observed=fseqs, expected=expected.freqs*length(simseqs) )
resids$z <- (resids$observed-resids$expected)/sqrt(resids$expected)
resids$p <- pbinom( resids$observed, size=length(simseqs), prob=expected.freqs )

cat("Beginning from ", as.character(initseq), ":\n")
print(resids)

if ( min( pbeta( min(resids$p), 1,nrow(resids) ) , pbeta( min(1-resids$p), 1,nrow(resids) ) ) < .01 ) {
    stop("Error: some high residuals there.")
}

