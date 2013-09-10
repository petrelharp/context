# took mtCDNApri.nuc from paml/examples/mtCDNA
#  and: removed extra spaces on first line (only one space between numbers)
#       made each sequence be on a single line
#       removed trailing text


require(Biostrings)

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

mtCDNA <- as( readDNAMultipleAlignment("mtCDNApri-sub.nuc", format="phylip" ), "DNAStringSet" )

lwin <- rwin <- 0
win <- 3
winlen <- lwin+win+rwin
boundary <- "none"; meanboundary <- 0
gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='')
load(gmfile)

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- lapply( list( hu.ch=c("human","chimpanzee"), ch.hu=c("chimpanzee","human"), go.hu=c("gorilla","human") ), function (x) {
        counttrans( rownames(projmatrix), colnames(projmatrix), mtCDNA[[x[1]]],  mtCDNA[[x[2]]],  lwin=lwin ) 
    } ) 
