#!/usr/bin/Rscript

library(contextual)
library(contextutils)
library(simcontext)

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
tlen <- .5
seqlen <- 1000
mutrates <- c(1,1)
selcoef <- c(-.5,.5)

initseq <- rinitseq(seqlen,bases)
simseqs <- simseq( seqlen, tlen, intiseq=initseq, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, bases=bases, fixfn=fixfn ) 

long <- 6
short <- 2

longpats <- getpatterns(long,bases=bases)
shortpats <- getpatterns(short,bases=bases)
long.counts <- counttrans( ipatterns=longpats, fpatterns=longpats, simseqs=simseqs, leftwin=0 )
short.counts <- counttrans( ipatterns=shortpats, fpatterns=shortpats, simseqs=simseqs, leftwin=0 )
proj.counts <- projectcounts( counts=long.counts, new.leftwin=0, new.longwin=short, new.shortwin=short, overlapping=TRUE )

all.counts <- data.frame( inpat=rownames(proj.counts@counts)[row(proj.counts@counts)], 
    outpat=colnames(proj.counts@counts)[col(proj.counts@counts)], 
    proj=as.vector(proj.counts@counts),
    stringsAsFactors=FALSE
    )
all.counts$short <- short.counts@counts[ cbind( match(all.counts$inpat,rownames(short.counts@counts)), match(all.counts$outpat,colnames(short.counts@counts)) ) ]

inseq <- as.character(simseqs$initseq)
inseq.vec <- substring( inseq, first=1:(seqlen-short+1), last=short:seqlen )
outseq <- as.character(simseqs$finalseq)
outseq.vec <- substring( outseq, first=1:(seqlen-short+1), last=short:seqlen )

all.counts$truth <- sapply( 1:nrow(all.counts), function (k) {
        sum( grepl( all.counts$inpat[k], inseq.vec ) & grepl( all.counts$outpat[k], outseq.vec ) )
    } )

stopifnot( all( all.counts$short == all.counts$truth ) )
stopifnot( sum(all.counts$proj)+long-short == sum(all.counts$short) )
stopifnot( all( ( all.counts$short - all.counts$proj ) <= long-short ) )

