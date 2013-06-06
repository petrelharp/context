#!/usr/bin/Rscript --vanilla
usage <- "\
Summarize results from mcmc runs.
"

if (!interactive()) {
    args <- commandArgs(TRUE)
    if (length(args)>=1) {
        simdir <- args[1]
    } else {
        simdir <- grep("sims$",list.dirs(recursive=FALSE),value=TRUE)
    }
}
if (!file.exists(simdir)) { stop(paste(simdir, "does not exist.\n")) }

simfiles <- list.files(simdir,"selsims-.*RData",full.names=TRUE)
basedirs <- gsub(".RData","",simfiles,fixed=TRUE)
mcmcdatafiles <- lapply(basedirs, function (x) list.files(path=x,pattern=".*-mcmc.*RData",full.names=TRUE) )
resultsfiles <- lapply(basedirs, function (x) list.files(path=x,pattern=".*-results.RData",full.names=TRUE) )
names(resultsfiles) <- names(simfiles) <- names(mcmcdatafiles) <- gsub(".*/","",basedirs)

# siminfo <- lapply( simfiles, function (x) { load(x); list(tlen=tlen/sum(tlen),muttime=mutrates*sum(tlen),selcoef=selcoef) } )
truths <- lapply( resultsfiles, function (x) { load(x); truth } )

mcmcruns <- lapply( mcmcdatafiles, function (dfs) {
        lapply(dfs, function (x) { tmp <- load(x);  list( lwin=lwin, win=win, rwin=rwin, mrun=mrun ) } ) } )
mcmcbatches <- lapply( mcmcruns, function (dfs) {
        do.call( rbind, lapply(dfs,function(x)x$mrun$batch) ) } )
hasmcmc <- (! sapply(mcmcbatches,is.null) )

mcmcinfo <- lapply( names(mcmcruns[hasmcmc]), function (x) { 
        bleh <- lapply( mcmcruns[[x]], function (y) {
                truth <- truths[[x]]
                # tlens <- siminfo[[x]]$tlen
                # muttimes <- siminfo[[x]]$muttime
                # selcoef <- siminfo[[x]]$selcoef
                z <- with(y, data.frame( 
                        fname=x, 
                        lwin=lwin, win=win, rwin=rwin, 
                        accept=mrun$accept, user.time=mrun$time[1], nbatch=mrun$nbatch, blen=mrun$blen, nspac=mrun$nspac
                        ) )
                ret <- with( y, do.call( cbind, c( list( z ), 
                        lapply( 1:ncol(mrun$batch), function (k) {
                            if (is.null(colnames(mrun$batch))) { colnames(mrun$batch) <- paste("param",1:ncol(mrun$batch),sep='') }
                            colnames(mrun$batch)[is.na(colnames(mrun$batch))] <- paste("param",(1:ncol(mrun$batch))[is.na(colnames(mrun$batch))], sep='')
                            thistruth <- truth[match(colnames(mrun$batch),names(truth))]
                            tmp <- data.frame( 
                                    q.truth=mean(mrun$batch[,k]<=thistruth[k]),
                                    q05=quantile(mrun$batch[,k],.05), q25=quantile(mrun$batch[,k],.25), 
                                    med=quantile(mrun$batch[,k],.50), 
                                    truth=thistruth[k],
                                    mean=mean(mrun$batch[,k]),
                                    q75=quantile(mrun$batch[,k],.75), q95=quantile(mrun$batch[,k],.95)
                                ) 
                            names(tmp) <- paste(colnames(mrun$batch)[k],names(tmp),sep="-")
                            return(tmp)
                        } )
                    ) ) )
                return(ret)
            } ) 
       do.call(rbind, lapply(bleh, function (x) { names(x) <- names(bleh[[1]]); x } ) )
   } )
mcmcinfo <- do.call( rbind, lapply( mcmcinfo, function (x) { names(x) <- names(mcmcinfo[[1]]); x } ) )
rownames(mcmcinfo) <- NULL

pdf(file="all-mcmc-runs.pdf",height=7,width=10)
# nmuts <- length(siminfo[[1]]$muttime)
# mutlabels <- paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) )
par(mar=c(4,4,1,0)+.1)
for (k in 1:ncol(mcmcbatches[[1]])) {
    layout(matrix(1:(2*(sum(hasmcmc)+1)%/%2),nrow=2))
    for (x in names(mcmcbatches[hasmcmc])) {
        truth <- truths[[x]][k]
        est <- hist( mcmcbatches[[x]][,k], breaks=50, plot=FALSE )
        xrange <- range(c(truth,est$breaks))
        xrange <- mean(xrange) + 1.2*(xrange-mean(xrange))
        plot( est, xlim=xrange, col=adjustcolor("blue",.4), xlab=x, border=grey(.6), main=if(k==1){names(truths[[x]])[k]} )
        abline(v=truth,lwd=2,col=adjustcolor('black',.75))
    }
}
dev.off()

write.table( format(mcmcinfo,digits=4), file=pipe("column -t"), row.names=FALSE, quote=FALSE )
