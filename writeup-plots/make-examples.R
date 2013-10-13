#/usr/bin/R --vanilla

source("../sim-context-fns.R")
source("../codon-inference-fns.R")


# EXAMPLE 1 in the talk
bases <- c("X","O")
mutpats <- list( 
    list( c("X","O") ), 
    list( c("O","X") ), 
    list( c("XO","OX") )
    ) 
mutrates <- c(1,1,1)
selpats <- list()
selcoef <- numeric(0)
fixfn <- function (...) { 1 }


initfreqs <- c(.2,.8)

seqlen <- 50
tlen <- .01
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
simseqs <- list( 
        simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases, count.trans=TRUE )
    )
for (k in 2:10) {
        simseqs <- c( simseqs, list(
                simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=simseqs[[k-1]]$finalseq, bases=bases, count.trans=TRUE )
                ) )
    }

lapply(simseqs,show.simseq,latex=TRUE,printit=TRUE)

# TASEP in the talk
bases <- c("X","O")
mutpats <- list( 
        list( c("XO","OX") )
    ) 
mutrates <- c(1)
selpats <- list()
selcoef <- numeric(0)
fixfn <- function (...) { 1 }

initfreqs <- c(.2,.8)

seqlen <- 30
tlen <- .5
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
simseqs <- list( 
        simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases, count.trans=TRUE )
    )
for (k in 2:10) {
        simseqs <- c( simseqs, list(
                simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=simseqs[[k-1]]$finalseq, bases=bases, count.trans=TRUE )
                ) )
    }


lapply(simseqs,show.simseq,latex=TRUE,printit=TRUE)


# ATGC in the talk
bases <- c("A","C","G","T")
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
mutrates <- c(rep(1,length(mutpats)-1),10)
selpats <- list()
selcoef <- numeric(0)
fixfn <- function (...) { 1 }

initfreqs <- c(.25,.25,.25,.25)

seqlen <- 30
tlen <- .05
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
simseqs <- list( 
        simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases, count.trans=TRUE )
    )


lapply(simseqs,show.simseq,latex=TRUE,printit=TRUE)
