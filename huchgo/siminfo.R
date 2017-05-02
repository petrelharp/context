#!/usr/bin/Rscript

usage <- "\
Print info about available simulations.\
Usage:\
   Rscript siminfo.R (directory with simulations)\
\
"

if (!interactive()) {
    args <- commandArgs(TRUE)
    simdir <- args[1]
}
if (!file.exists(simdir)) { stop(paste(simdir, "does not exist.\n")) }


# available simulated sequences
simfiles <- list.files(simdir,"*.RData",full.names=TRUE)
siminfo <- lapply(simfiles, function (x) {
        basedir <- gsub(".RData","",x,fixed=TRUE)
        tmp <- load(x)
        y <- list(date=now,seqlen=seqlen, tlen=tlen, file=x) 
        pnames <- names(y)
        if ("mutrates"%in%tmp) { 
            y <- c( y, as.list(mutrates) )
            pnames <- c( pnames, paste("mutrate",seq_along(mutrates),sep='') )
        }
        if ("selcoef"%in%tmp & length(selcoef)>0) { 
            y <- c( y, as.list(selcoef) )
            pnames <- c( pnames, paste("selcoef",seq_along(selcoef),sep='') )
        }
        y <- c( y, list( stringsAsFactors=FALSE ) )
        y <- do.call( data.frame, y )
        colnames(y) <- pnames
        estimfiles <- list.files(basedir,"*results.tsv",full.names=TRUE)
        for (ef in estimfiles) {
            wins <- gsub("-results.tsv","",gsub(".*win-","",ef))
            estim <- read.table(ef,sep="\t")
            estim <- estim[match(c("ans","random.ans"),rownames(estim))[1],setdiff(colnames(estim),"likfun"),drop=FALSE]
            colnames(estim) <- paste( wins, colnames(estim), sep=".")
            y <- cbind(y, estim)
        }
        return(y)
    } )
allnames <- unique( as.vector( unlist( sapply( siminfo, colnames ) ) ) )
siminfo <- lapply( siminfo, function (x) { for (y in setdiff(allnames, colnames(x))) { x[[y]] <- NA }; x[allnames] } )
siminfo <- do.call(rbind, siminfo)
siminfo$meandist <- siminfo$tlen * rowSums(siminfo[,grep("mutrate",colnames(siminfo)),drop=FALSE])
siminfo <- siminfo[ order(siminfo$meandist), ]

# switch to muttime
mutvars <- grep("mutrate",colnames(siminfo),value=TRUE)
for (mm in mutvars) {
    tt <- gsub("rate","time",mm)
    siminfo[[tt]] <- siminfo$tlen * siminfo[[mm]]
}
mtvars <- grep("muttime",colnames(siminfo),value=TRUE)
mtvars <- mtvars[order(gsub(".*\\.","",mtvars))]
selvars <- grep("selcoef",colnames(siminfo),value=TRUE)
selvars <- selvars[order(gsub(".*\\.","",selvars))]
siminfo <- siminfo[, c( setdiff(colnames(siminfo),c(mtvars,selvars,mutvars)), mtvars, selvars ) ]

write.table(format(siminfo,digits=3),file=pipe("column -t -s '\t'"),quote=FALSE,row.names=FALSE,sep="\t")

