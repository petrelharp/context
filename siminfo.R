#!/usr/bin/Rscript --vanilla

usage <- "\
Print info about available simulations.\
Usage:\
   Rscript siminfo.R (directory with simulations)\
\
"

args <- commandArgs(TRUE)
if (length(args)<1) {
    simdir <- "sims/"
} else {
    simdir <- args[1]
}


# available simulated sequences
simfiles <- list.files(simdir,"*.RData",full.names=TRUE)
siminfo <- do.call(rbind, lapply(simfiles, function (x) {
        tmp <- load(x)
        y <- list(date=now,seqlen=seqlen, tlen=tlen, file=x) 
        pnames <- names(y)
        if ("mutrates"%in%tmp) { 
            y <- c( y, as.list(mutrates) )
            pnames <- c( pnames, paste("mutrate",seq_along(mutrates),sep='') )
        }
        if ("selcoef"%in%tmp) { 
            y <- c( y, as.list(selcoef) )
            pnames <- c( pnames, paste("selcoef",seq_along(selcoef),sep='') )
        }
        y <- c( y, list( stringsAsFactors=FALSE ) )
        y <- do.call( data.frame, y )
        colnames(y) <- pnames
        y
    } ) )
siminfo$meandist <- siminfo$tlen * colSums(siminfo[,grep("mutrate",colnames(siminfo)),drop=FALSE])
siminfo <- siminfo[ order(siminfo$meandist), ]

write.table(as.matrix(siminfo),file=pipe("column -t -s '\t'"),quote=FALSE,row.names=FALSE,sep="\t")

