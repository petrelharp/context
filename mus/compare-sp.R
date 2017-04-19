#!/usr/bin/Rscript

library(contextual)
library(contextutils)

aligns <- c( "mm9ornAna1", "mm9hg19", "mm9galGal3", "mm9oryCun2", "mm9rn5" )

shortwin <- 3
leftwin <- rightwin <- 1
longwin <- shortwin+leftwin+rightwin
bases <- c("A","T","C","G")
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=list(), boundary='none' )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )

count.files <- outer( aligns, c("5.1.counts","rev.5.1.counts"), paste, sep='/' )
names(count.files) <- gsub("\\/.*","",count.files)
count.files <- count.files[order(names(count.files))]
counts <- lapply(count.files, function (ifile) {
        count.table <- read.table(ifile,header=TRUE,stringsAsFactors=FALSE)
        counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
        rownames(counts) <- rownames(genmatrix)
        colnames(counts) <- colnames(projmatrix)
        stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
        counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
        return(counts)
    } )
simple.counts <- lapply(counts, countmuts, mutpats=mutpats, leftwin=leftwin)
simple.table <- data.frame( do.call(rbind, lapply(simple.counts,t) ) )
names(simple.table) <- c("numerator","denominator")
simple.table$mut <- factor( sapply(simple.counts,colnames) )
simple.table$sp <- factor( names(simple.counts)[rep(seq_along(simple.counts),sapply(simple.counts,ncol))] )

divergences <- sapply( counts, divergence, leftwin=leftwin )

with(simple.table, { plot( denominator, numerator, col=sp, pch=as.numeric(mut) );
        legend("topright", col=seq_along(levels(sp)), legend=levels(sp), pch=1);
        legend("topleft", legend=levels(mut), pch=seq_along(levels(mut)) )
    } )

## estimates
mle.tables <- lapply( sapply(paste(aligns,"/3.1.counts-results",sep=''),list.files,pattern="*.tsv",full.names=TRUE),
    read.csv, header=TRUE, sep='\t' )
names(mle.tables) <- aligns
mles <- data.frame( t( sapply(mle.tables,function (x) as.matrix(x)["mle",]) ) )
colnames(mles) <- gsub("..","->",gsub("...","|",gsub("tmut..","",colnames(mles)),fixed=TRUE),fixed=TRUE)

mles <- mles[order(mles$branchlen),]
matplot(mles[,1:18],type='l')
legend('topleft',legend=colnames(mles)[1:18],col=1:6,lty=1:5)

