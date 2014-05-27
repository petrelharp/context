# make colored pairplots of combined mcmc runs

args <- commandArgs(TRUE)
simdir <- args[1]

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))
require(mcmc)

if (interactive()) {
    simdir <- "cpg-tree-sims"
}

# infile <- "cpg-tree-sims/selsims-2013-06-03-13-17-0356374.RData"
for (infile in list.files(simdir,"*RData",full.names=TRUE)) {

    basedir <- gsub(".RData","",infile,fixed=TRUE)
    load(infile)

    lwin <- rwin <- 2; win <- 1
    basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
    datafile <- paste( basename ,"-results.RData",sep='')
    plotfile <- paste( basename ,"-plot-all-mcmcs",sep='')
    mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
    mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

    load(datafile)

    varnames <- switch( simdir,
            "cpg-tree-sims" = c(
                    "branchlen", 
                    paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), 
                    names(initfreqs)[-length(initfreqs)] 
                ),
            "cpg-sims" = c(paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ) ),
            "ising-sims" = c("lambda","beta","gamma")
        )


    all.mrun <- do.call( rbind, lapply( seq_along(mcmcdatafiles), function (k) {
                load(mcmcdatafiles[[k]])
                tmp <- data.frame( mrun$batch )
                names(tmp) <- varnames
                tmp$mcmcnum <- k
                tmp
            } ) )

    subseq <- seq(1,nrow(all.mrun),by=max(1,floor(nrow(all.mrun)/1000)))

    x <- rbind( all.mrun[subseq,(-1)*match("mcmcnum",colnames(all.mrun))], truth ) 

    # # color by run
    # cols <- c( adjustcolor(rainbow(7),.2)[all.mrun$mcmcnum[subseq]], "black" )
    # color by time
    pdf(file=paste(plotfile,".pdf",sep=''),width=12,height=12,pointsize=10)
    cols <- c( adjustcolor(rainbow(64),.2)[1+floor(65*(1:(nrow(x)-1))/nrow(x))], adjustcolor("black",.5) )
    pairs( x, col=cols, pch=20, cex=c(rep(.5,nrow(x)-1),2) )
    dev.off()
}
