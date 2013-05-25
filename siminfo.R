#!/usr/bin/Rscript --vanilla

usage <- "\
Print info about available simulations.\
Usage:\
   Rscript siminfo.R (directory with simulations)\
\
"

args <- commandArgs(TRUE)
if (length(args)<1) {
    simdir <- "ising-sims/"
} else {
    simdir <- args[1]
}


# available simulated sequences
simdir <- "ising-sims"
simfiles <- list.files(simdir,"*.RData",full.names=TRUE)
siminfo <- do.call(rbind, lapply(simfiles, function (x) {
        load(x)
        y <- do.call( data.frame, c(list(date=now,seqlen=seqlen, tlen=tlen, file=x), as.list(mutrates), as.list(selcoef), stringsAsFactors=FALSE) )
        colnames(y) <- c( colnames(y)[1:4], paste("mutrate",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
        y
    } ) )
siminfo$meandist <- siminfo$tlen * colSums(siminfo[,grep("mutrate",colnames(siminfo)),drop=FALSE])
siminfo <- siminfo[ order(siminfo$meandist), ]

write.table(as.matrix(siminfo),file=pipe("column -t -s '\t'"),quote=FALSE,row.names=FALSE,sep="\t")

