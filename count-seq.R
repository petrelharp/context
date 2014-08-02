
#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Count paired tuples from the Rdata file containing simulated (or other) sequence,\
of the following form: \
    WWWWWWWWW \
    llMMMMrrr \
where the  length of the W's is 'winlen', the length of the l's is 'lwin', and the length of the 'M's is 'win'.
 \
If --revcounts, count in the other direction, i.e. \
    llMMMMrrr \
    WWWWWWWWW \
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", help=".RData file containing simulation." ),
        make_option( c("-o","--outfile"), type="character", help="File to direct output to."),
        make_option( c("-w","--winlen"), type="integer", help="Size of long window." ),
        make_option( c("-s","--win"), type="integer", help="Size of short window." ),
        make_option( c("-l","--lwin"), type="integer", help="Size of offset of short window from the left."),
        make_option( c("-r","--revcounts"), action="store_true", default="FALSE", help="Count reversed?")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$outfile)) { outfile <- paste( gsub(".RData","",opt$infile), if(opt$revcounts){"-rev"}else{""}, ".counts", sep="") }
attach(opt)

source("../context-inference-fns.R")
source("../sim-context-fns.R")

load(infile)

longpats <- getpatterns(winlen,bases)
shortpats <- getpatterns(win,bases)

# this returns a matrix
counts <- if (revcounts) {
    counttrans( longpats, shortpats, simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin )
} else {
    counttrans( longpats, shortpats, simseqs[[1]]$finalseq, simseqs[[1]]$initseq, lwin=lwin )
}

countframe <- data.frame( reference=rownames(counts)[row(counts)], 
                         derived=colnames(counts)[col(counts)],
                         count=as.vector(counts)
                         )

write.table(countframe, file=outfile, row.names=FALSE, sep=' ', quote=FALSE)
