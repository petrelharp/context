#!/usr/bin/Rscript --vanilla
require(optparse)
require(colorspace)

usage <- "\
Plot simple statistics for count files in some directory.
"

option_list <- list(
    # input/output
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-o","--outfile"), type="character", default="", help="Direct output to this file. [default \"simple-stats-shortwin-leftwin-rightwin.pdf\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
    # context and pattern size
        make_option( c("-w","--shortwin"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lrwin"), type="integer", default=1, help="Size of left-hand and right-hand context. [default \"%default\"]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$basedir)) { stop("No input directory  Run\n  bcells-inference.R -h\n for help.\n") }
attach(opt)

leftwin <- rightwin <- lrwin
longwin <- shortwin + leftwin + rightwin

if (opt$outfile=="") { outfile <- with(opt, paste(basedir,"/simple-stats",shortwin,"-",leftwin,"-",rightwin,".pdf",sep='')) }

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))

countfiles <- list.files(basedir,paste("*",longwin,lrwin,"counts",sep="."), full.names=TRUE)
names(countfiles) <- gsub(".*chr([0-9XY]*).*","chr\\1",countfiles)
countfiles <- countfiles[ order( as.numeric( gsub("chr","",names(countfiles)) ) ) ]
if (length(countfiles)==0) { stop("No input files.") }

boundary <- "none"
meanboundary <- 0
patlen <- 1
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,patlen,sep="-"),".RData",sep='') }
if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop("Can't find generator matrix in ", gmfile, " -- provide file name exactly?")
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )

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
        adhoc <- countmuts(counts=counts,mutpats=mutpats,lightwin=lightwin)
        adhoc <- adhoc[1,]/adhoc[2,]
    } )
rownames(adhoc) <- mutnames <- sapply(mutpats,function(x){paste(sapply(x,paste,collapse='->'),collapse='/')})

mean.counts <- Reduce("+",counts)/length(counts)

pdf(file=outfile, width=10, height=7.5, pointsize=10)
layout(matrix(1:min(ncol(counts),6),nrow=2))
par(mar=c(5,2,1,1)+.1)
cols <- rainbow_hcl(24)
matplot( adhoc, type='l', ylab='mutation rate', xlab='', xaxt='n', col=cols )
axis(1, at=1:nrow(adhoc), labels=mutnames, las=3 )
for (k in 1:ncol(mean.counts)) {
    lord <- order( mean.counts[,k] )
    plot( mean.counts[lord,k], xaxt='n', xlab='', main=colnames(mean.counts)[k], log='y', ylim=1+range(unlist(sapply(counts,as.numeric))), ylab='' )
    axis(1,at=1:nrow(mean.counts),labels=rownames(mean.counts)[lord],las=3)
    for (j in seq_along(counts)) {
                points( counts[[j]][lord,k], col=adjustcolor(cols[j],.5), cex=.5, pch=20 ) 
                lines( lowess(seq_along(lord),counts[[j]][lord,k],f=.05), col=adjustcolor(cols[j],.5), lwd=2 )
            }
    legend("topleft",legend=names(counts),lty=1,col=cols)
}
dev.off()
