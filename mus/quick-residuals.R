#!/usr/bin/R

source("../context-inference-fns.R")

long.win <- list( lwin=2, rwin=2, win=3, winlen=7 )
# long.win <- list( lwin=1, rwin=1, win=5, winlen=7 )

bases <- c("A","T","C","G")
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
selpats <- list()
fixfn <- function (...) { 1 }
boundary <- 'none'
meanboundary <- 0

long.win$genmatrix <- with(long.win, makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=list(), boundary=boundary ) )
long.win$projmatrix <- with(long.win, collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin ) )

# all kmer -> kmer changes
long.mutpats <- lapply(1:3, function (k) {
    kmers <- getpatterns(k)
    c( apply(combn(kmers,2),2,list),  
            apply(combn(kmers,2)[2:1,],2,list)
        )
} )


# gmfile <- paste(paste("genmatrices/genmatrix",commandArgs[1],'none',0,sep="-"),".RData",sep='') 
# load(gmfile)
# long.genmat <- genmatrix

setwd("mm9rn5")

getcounts <- function (winlen,win,genmatrix,projmatrix) {
    infile <- paste(winlen,win,"counts",sep='.')
    revfile <- paste("rev.",infile,sep='')
    counts <- lapply( list(infile,revfile), function (ifile) {
            count.table <- read.table(ifile,header=TRUE,stringsAsFactors=FALSE)
            counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
            rownames(counts) <- rownames(genmatrix)
            colnames(counts) <- colnames(projmatrix)
            stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
            counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
            return(counts)
        } )
    return(counts)
}

long.win$counts <- with(long.win, getcounts(winlen,win,genmatrix,projmatrix))


load("3.1.counts-results/win-3-1-1-results.RData")
load("3.1.counts-results/win-1-1-1-mcmc-2.RData")

if (mrun$nbatch * mrun$blen > 1e4) {
    estimates <- rbind(estimates, 
           final=c(mrun$final,1-sum(mrun$final[length(mrun$final)-0:2]),NA),
           mean=c(colMeans(mrun$batch),1-sum(colMeans(mrun$batch[,length(mrun$final)-0:2])),NA))
}

expected <- function (x,win,lwin,rwin,genmatrix,projmatrix,counts) {
            x <- as.numeric(x)
            branchlens <- c(x[1],1-x[1])
            mutrates <- x[1+(1:nmuts)]
            initfreqs <- x[1+nmuts+(1:nfreqs)]
            initfreqs <- initfreqs/sum(initfreqs)
            names(initfreqs) <- bases
            list( 
                    predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[1]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), 
                            genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=rev(branchlens) ),
                    predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[2]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), 
                            genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=branchlens )
                )
}

all.expected <- with( long.win, apply(estimates,1,expected,win,lwin,rwin,genmatrix,projmatrix,counts) )
names(all.expected) <- rownames(estimates)

# all.resids <- with(long.win, lapply( all.expected, function (x) mapply(function(u,v) (u-v)/sqrt(u),x,counts) ) )
all.resids <- with(long.win, lapply( all.expected, function (x) mapply(function(u,v) qnorm(ppois(as.numeric(v),as.numeric(u))), x, counts, SIMPLIFY=FALSE)) )
for (i in seq_along(all.resids)) for (j in seq_along(all.resids[[i]])) { { with(long.win, { dim(all.resids[[i]][[j]]) <- dim(counts[[1]]); dimnames(all.resids[[i]][[j]]) <- dimnames(counts[[1]]) } ) } }

resid.counts <- with(long.win, lapply( all.resids[['mle']], countmuts, mutpats=long.mutpats[[2]], lwin=lwin) )

resid.counts <- with(long.win, lapply( long.mutpats, function (mutpats) { lapply( all.resids[['mle']], countmuts, mutpats=mutpats, lwin=lwin) } ) )

save(list=ls(),file=with(long.win, paste("resids-",win,"-",lwin,".RData",sep='')))


layout(matrix(1:6,nrow=2))
for (j in 2:4) {
    with(long.win, {
        plot( as.numeric( counts[[1]]), as.numeric(all.resids[[j]][[1]]) );
        plot( as.numeric( counts[[2]]), as.numeric(all.resids[[j]][[2]]) )
    } )
}


layout(matrix(1:6,nrow=3))
lapply( lapply(all.resids[2:4],"[[",1), matplot, type='l', col=adjustcolor(1:6,.2) )
lapply( lapply(all.resids[2:4],"[[",2), matplot, type='l', col=adjustcolor(1:6,.2) )


###

cwin <- min(2,win); lrcwin <- min(1,lwin,rwin)
subcounts <- lapply( counts, function (x) 
        projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )
all.subexpected <- lapply( all.expected, lapply, function (x)
        projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )
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
