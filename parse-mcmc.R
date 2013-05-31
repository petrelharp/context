#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Summarize results from mcmc runs.
"

option_list <- list(
        make_option( c("-i","--simdir"), type="character", default=grep("sims$",list.dirs(recursive=FALSE),value=TRUE), help="directory containing " ),
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

simfiles <- list.files(simdir,"selsims-.*RData",full.names=TRUE)
basedirs <- gsub(".RData","",simfiles,fixed=TRUE)
mcmcdatafiles <- lapply(basedirs, function (x) list.files(path=x,pattern=".*-mcmc.*RData",full.names=TRUE) )
names(simfiles) <- names(mcmcdatafiles) <- gsub(".*/","",basedirs)

siminfo <- lapply( simfiles, function (x) { load(x); list(muttime=mutrates*tlen) } )

mcmcruns <- lapply( mcmcdatafiles, function (dfs) {
        lapply(dfs, function (x) { tmp <- load(x);  list( lwin=lwin, win=win, rwin=rwin, mrun=mrun ) } ) } )
mcmcbatches <- lapply( mcmcruns, function (dfs) {
        do.call( c, lapply(dfs,function(x)x$mrun$batch) ) } )

mcmcinfo <- do.call(rbind, lapply( names(mcmcruns), function (x) { do.call( rbind, lapply( mcmcruns[[x]], function (y) {
                z <- with(y, data.frame( 
                        fname=x, 
                        lwin=lwin, win=win, rwin=rwin, 
                        accept=mrun$accept, user.time=mrun$time[1], nbatch=mrun$nbatch, blen=mrun$blen, nspac=mrun$nspac,
                        q.truth=mean(mrun$batch<=siminfo[[x]]$muttime),
                        q05=quantile(mrun$batch,.05), q25=quantile(mrun$batch,.25), 
                        med=quantile(mrun$batch,.50), 
                        truth=siminfo[[x]]$muttime,
                        mean=mean(mrun$batch),
                        q75=quantile(mrun$batch,.75), q95=quantile(mrun$batch,.95)
                        ) )
            } ) ) } ) )
rownames(mcmcinfo) <- NULL

pdf(file="all-mcmc-runs.pdf",height=7,width=10)
layout(matrix(seq_along(simfiles),nrow=2))
par(mar=c(4,4,0,0)+.1)
for (x in names(mcmcbatches)) {
    truth <- siminfo[[x]]$muttime
    est <- hist( mcmcbatches[[x]], breaks=50, plot=FALSE )
    xrange <- range(c(truth,est$breaks))
    xrange <- mean(xrange) + 1.2*(xrange-mean(xrange))
    plot( est, xlim=xrange, col=adjustcolor("blue",.4), main='', xlab=x )
    abline(v=truth,lwd=2,col=adjustcolor('black',.75))
}
dev.off()

print( mcmcinfo, file=pipe("column -t") )
