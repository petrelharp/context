#!/usr/bin/Rscript
# Test inference from simulation

usage <- "\
Usage:\
   Rscript sim-seq.R seqlen tlen\
"

args <- commandArgs(TRUE)
if (length(args)<2) {
    stop(usage)
} else {
    seqlen <- as.numeric(args[1])  # e.g. 1e4
    tlen <- as.numeric(args[2]) # total length 6e7 gives lots of transitions (but still signal); 1e7 not so many
}

source("codons.R")
source("sim-context-fns.R")
source("codon-inference-fns.R")

# maximum size of pattern (for simulation)
patlen <- 2
mutpats <- list( 
        fwd=combn( bases, 2, simplify=FALSE ), 
        rev=lapply( combn( bases, 2, simplify=FALSE ), rev),
        gc=list( c("GC","TC") ) ,
        ga=list( c("GA","CA") ) ,
        ta=list( c("AT","GT") ) 
    ) 
mutrates <- runif( length(mutpats) )*1e-8
selpats <- c( "[GC]", "[AT]" )
selcoef <- runif( length(selpats) )*1e-3
# other params
Ne <- c(1e4,1e4)
branchlens <- c(1,1)

initfreqs <- c(.3,.3,.2,.2)
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- lapply(1:2, function (k) simseq( seqlen, tlen*branchlens[k], patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, Ne=Ne[k], initseq=initseq ) ) 
    )

thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- format(Sys.time(), "%Y-%m-%d-%H-%M")
save( thisone, now, patlen, mutpats, selpats, selcoef, Ne, seqlen, tlen, branchlens, initfreqs, simseqs, file=paste(now,thisone,"sims/selsims.RData",sep='') )


