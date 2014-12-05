#!/usr/bin/Rscript
# usage:
#  quick-residuals.R (length of window) (files to load)

gmfile <- paste(paste("genmatrices/genmatrix",commandArgs[1],'none',0,sep="-"),".RData",sep='') 
load(gmfile)
big.genmat <- genmatrix

for (x in commandArgs(TRUE)[-1]) { load(x) }

if (mrun$nbatch * mrun$blen > 1e4) {
    estimates <- rbind(estimates, 
           final=c(mrun$final,1-sum(mrun$final[length(mrun$final)-0:2]),NA),
           mean=c(colMeans(mrun$batch),1-sum(colMeans(mrun$batch[,length(mrun$final)-0:2])),NA))
}

counts <- lapply( list(infile,revfile), function (ifile) {
        count.table <- read.table(ifile,header=TRUE,stringsAsFactors=FALSE)
        counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
        rownames(counts) <- rownames(genmatrix)
        colnames(counts) <- colnames(projmatrix)
        stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
        counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
        return(counts)
    } )

expected <- function (x) {
            x <- as.numeric(x)
            branchlens <- c(x[1],1-x[1])
            mutrates <- x[1+(1:nmuts)]
            initfreqs <- x[1+nmuts+(1:nfreqs)]
            initfreqs <- initfreqs/sum(initfreqs)
            names(initfreqs) <- bases
            list( 
                    predicttreecounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[1]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), 
                            genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=rev(branchlens) ),
                    predicttreecounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[2]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), 
                            genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=branchlens )
                )
}

all.expected <- apply(estimates,1,expected)
names(all.expected) <- rownames(estimates)


cwin <- min(2,shortwin); lrcwin <- min(1,leftwin,rightwin)
subcounts <- lapply( counts, function (x) 
        projectcounts( x, new.shortwin=cwin, new.leftwin=lrcwin, new.longwin=2*lrcwin+cwin, counts=x ) )
all.subexpected <- lapply( all.expected, lapply, function (x)
        projectcounts( x, new.shortwin=cwin, new.leftwin=lrcwin, new.longwin=2*lrcwin+cwin, counts=x ) )
all.subresids <- lapply( all.subexpected, function (x) mapply(function(u,v) (u-v)/sqrt(v),x,subcounts) )

layout(matrix(seq_along(subcounts)))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
for (j in seq_along(counts)) {
    z <- sapply( lapply( all.subresids, "[[", j ), as.vector )
    rownames(z) <- paste( rownames(subcounts[[j]])[row(subcounts[[j]])], colnames(subcounts[[j]])[col(subcounts[[j]])], sep="->" )
    plot( 0, type='n', xlab="", ylab="normalized residuals", xlim=c(0,ncol(z)+1), ylim=range(z), xaxt='n' )
    text( jitter(col(z),fac=2), z, labels=rownames(z), col=col(z) )
    axis(1, at=1:ncol(z), labels=colnames(z) )
}
