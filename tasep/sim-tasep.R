#!/usr/bin/Rscript
# Simulate up the process

usage <- "\
Simulate from the process.\
\
Usage:\
   Rscript sim-tasep.R seqlen tlen [Xfreq] \
"

args <- commandArgs(TRUE)
if (length(args)<2) {
    stop(usage)
} else {
    seqlen <- as.numeric(args[1])  # e.g. 1e4
    tlen <- as.numeric(args[2]) # total length 6e7 gives lots of transitions (but still signal); 1e7 not so many
    if (length(args)>=3) {
        Xfreq <- as.numeric(args[3]) # density of Xs
    } else {
        Xfreq <- 0.5
    }
}

simdir <- "tasep-sims/"

bases <- c("X","O")

source("../sim-context-fns.R")
source("../codon-inference-fns.R")

# maximum size of pattern (for simulation)
patlen <- 2
mutpats <- list( 
    list( c("XO","OX") )
    ) 
mutrates <- 2 * runif( length(mutpats) )

initfreqs <- c(Xfreq,1-Xfreq)
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- list(
                simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, mutrates=mutrates, selpats=c(), selcoef=numeric(0), Ne=1, initseq=initseq, bases=bases )
            )
    )

thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- Sys.time()
save( thisone, now, patlen, mutpats, mutrates, seqlen, tlen, initfreqs, simseqs, file=paste(simdir,format(now,"%Y-%m-%d-%H-%M"),thisone,"selsims.RData",sep='') )


