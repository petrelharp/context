#!/usr/bin/Rscript
# Simulate up the process

usage <- "\
Simulate from the process.\
\
Usage:\
   Rscript sim-ising.R seqlen tlen interaction external_field [Xfreq] \
"

args <- commandArgs(TRUE)
if (length(args)<2) {
    stop(usage)
} else {
    seqlen <- as.numeric(args[1])  # e.g. 1e4
    tlen <- as.numeric(args[2]) # total length 6e7 gives lots of transitions (but still signal); 1e7 not so many
    interaction <- as.numeric(args[3])
    external.field <- as.numeric(args[4])
    if (length(args)>=5) {
        Xfreq <- as.numeric(args[5]) # density of Xs
    } else {
        Xfreq <- 0.5
    }
}

# identifiers
thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- Sys.time()

simdir <- "ising-sims/"
if (!file.exists(simdir)) { dir.create(simdir,recursive=TRUE) }

bases <- c("X","O")

source("../sim-context-fns.R")
source("../codon-inference-fns.R")

# maximum size of pattern (for simulation)
patlen <- 2
mutpats <- list( 
    list( c("O","X"), c("X","O") )
    ) 
mutrates <- c(1)
selpats <- list(
        c("OX","XO"),
        c("X")
    )
selcoef <- c(interaction,external.field)

fixfn <- ising.fixfn

initfreqs <- c(Xfreq,1-Xfreq)
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- list(
                simseq( seqlen=seqlen, tlen=tlen, patlen=patlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases )
            )
    )

save( thisone, now, bases, patlen, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, initfreqs, simseqs, file=paste(simdir,"selsims-",format(now,"%Y-%m-%d-%H-%M"),"-",thisone,".RData",sep='') )


