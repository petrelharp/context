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

mcmcinfo <- do.call(rbind, lapply( names(mcmcruns[hasmcmc]), function (x) { do.call( rbind, lapply( mcmcruns[[x]], function (y) {
                truth <- truths[[x]]
                # tlens <- siminfo[[x]]$tlen
                # muttimes <- siminfo[[x]]$muttime
                # selcoef <- siminfo[[x]]$selcoef
                z <- with(y, data.frame( 
                        fname=x, 
                        lwin=lwin, win=win, rwin=rwin, 
                        accept=mrun$accept, user.time=mrun$time[1], nbatch=mrun$nbatch, blen=mrun$blen, nspac=mrun$nspac
                        ) )
                with( y, do.call( cbind, c( list( z ), 
                        lapply( seq_along(truth), function (k) {
                            thistruth <- truth[match(colnames(mrun$batch),names(truth))]
                            data.frame( q.truth=mean(mrun$batch[,k]<=thistruth[k]),
                            q05=quantile(mrun$batch[,k],.05), q25=quantile(mrun$batch[,k],.25), 
                            med=quantile(mrun$batch[,k],.50), 
                            truth=thistruth[k],
                            mean=mean(mrun$batch[,k]),
                            q75=quantile(mrun$batch[,k],.75), q95=quantile(mrun$batch[,k],.95)
                        ) } )
                    ) ) )
            } ) ) } ) )
rownames(mcmcinfo) <- NULL

pdf(file="all-mcmc-runs.pdf",height=7,width=10)
nmuts <- length(siminfo[[1]]$muttime)
mutlabels <- paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) )
par(mar=c(4,4,1,0)+.1)
for (k in seq_along(siminfo[[1]]$muttime)) {
    layout(matrix(1:(2*(sum(hasmcmc)+1)%/%2),nrow=2))
    for (x in names(mcmcbatches[hasmcmc])) {
        truth <- siminfo[[x]]$muttime[k]
        est <- hist( mcmcbatches[[x]][,k], breaks=50, plot=FALSE )
        xrange <- range(c(truth,est$breaks))
        xrange <- mean(xrange) + 1.2*(xrange-mean(xrange))
        plot( est, xlim=xrange, col=adjustcolor("blue",.4), xlab=x, border=grey(.6), main=if(k==1){mutlabels[k]} )
        abline(v=truth,lwd=2,col=adjustcolor('black',.75))
    }
}
sellabels <- paste("sel:", unlist(sapply(selpats,paste,collapse=" | ")) )
for (k in seq_along(siminfo[[1]]$selcoef)) {
    layout(matrix(1:(2*(sum(hasmcmc)+1)%/%2),nrow=2))
    for (x in names(mcmcbatches[hasmcmc])) {
        truth <- siminfo[[x]]$selcoef[k]
        est <- hist( mcmcbatches[[x]][,k+nmuts], breaks=50, plot=FALSE )
        xrange <- range(c(truth,est$breaks))
        xrange <- mean(xrange) + 1.2*(xrange-mean(xrange))
        plot( est, xlim=xrange, col=adjustcolor("red",.4), xlab=x, border=grey(.6), main=if(k==1){sellabels[k]} )
        abline(v=truth,lwd=2,col=adjustcolor('black',.75))
    }
}
dev.off()

write.table( format(mcmcinfo,digits=4), file=pipe("column -t"), row.names=FALSE, quote=FALSE )
