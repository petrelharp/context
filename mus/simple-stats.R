#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Plot simple statistics for count files in some directory.
"

option_list <- list(
    # input/output
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-o","--outfile"), type="character", default="", help="Direct output to this file. [default \"simple-stats-win-lwin-rwin.pdf\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
    # context and pattern size
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lrwin"), type="integer", default=2, help="Size of left-hand and right-hand context. [default \"%default\"]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$basedir)) { stop("No input directory  Run\n  bcells-inference.R -h\n for help.\n") }
attach(opt)

lwin <- rwin <- lrwin
winlen <- win - lwin - rwin

if (opt$outfile=="") { outfile <- with(opt, paste(basedir,"simple-stats-",win,"-",lwin,"-",rwin,".pdf",sep='')) }

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))

boundary <- "none"
meanboundary <- 0
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,patlen,sep="-"),".RData",sep='') }
if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop("Can't find generator matrix in ", gmfile, " -- provide file name exactly?")
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, time="gamma" )

countfiles <- list.files(basedir,paste("*",win,lrwin,"counts",sep="."), full.names=TRUE)
names(countfiles) <- gsub(".*chr([0-9]*).*","chr\\1",countfiles)
countfiles <- countfiles[ order( as.numeric( gsub("chr","",names(countfiles)) ) ) ]
counts <- lapply(countfiles, function (x) {
        count.table <- read.table(x,header=TRUE,stringsAsFactors=FALSE)
        counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
        rownames(counts) <- rownames(genmatrix)
        colnames(counts) <- colnames(projmatrix)
        stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
        counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
        return( counts )
    } )

adhoc <- sapply( counts, function (counts) {
        adhoc <- countmuts(counts=counts,mutpats=mutpats,lwin=lwin)
        adhoc <- adhoc[1,]/adhoc[2,]
    } )
rownames(adhoc) <- mutnames <- sapply(mutpats,function(x){paste(sapply(x,paste,collapse='->'),collapse='/')})

pdf(file=outfile, width=10, height=7.5, pointsize=10)
layout(matrix(1:min(ncol(counts),6),nrow=2))
par(mar=c(5,2,1,1)+.1)
matplot( adhoc, type='l', ylab='mutation rate', xlab='', xaxt='n' )
axis(1, at=1:nrow(adhoc), labels=mutnames, las=3 )
for (k in 1:ncol(counts)) {
    lord <- order( subexpected[,k] )
    plot( counts[lord,k], xaxt='n', xlab='', main=colnames(counts)[k], log='y', ylim=1+range(c(as.matrix(counts),as.matrix(subexpected))), ylab='' )
    axis(1,at=1:nrow(counts),labels=rownames(counts)[lord],las=3)
    lines(subexpected[lord,k],col='red')
    legend("topleft",legend='fitted',lty=1,col='red')
}
dev.off()
