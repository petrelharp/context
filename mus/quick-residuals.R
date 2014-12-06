#!/usr/bin/R

source("../context-inference-fns.R")

# long.win <- list( leftwin=2, rightwin=2, shortwin=3, longwin=7 )
long.win <- list( leftwin=1, rightwin=1, shortwin=5, longwin=7 )

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

long.win$genmatrix <- with(long.win, makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=list(), boundary=boundary ) )
long.win$projmatrix <- with(long.win, collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin ) )

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

getcounts <- function (longwin,shortwin,genmatrix,projmatrix) {
    infile <- paste(longwin,shortwin,"counts",sep='.')
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

long.win$counts <- with(long.win, getcounts(longwin,shortwin,genmatrix,projmatrix))


load("3.1.counts-results/win-3-1-1-results.RData")
load("3.1.counts-results/win-1-1-1-mcmc-2.RData")
# load("3.1.counts-dual-results/win-1-1-1-results.RData")

if (mrun$nbatch * mrun$blen > 1e4) {
    estimates <- rbind(estimates, 
           final=c(mrun$final,1-sum(mrun$final[length(mrun$final)-0:2]),NA),
           mean=c(colMeans(mrun$batch),1-sum(colMeans(mrun$batch[,length(mrun$final)-0:2])),NA))
}

expected <- function (x,shortwin,leftwin,rightwin,genmatrix,projmatrix,counts) {
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

all.expected <- with( long.win, apply(estimates,1,expected,shortwin,leftwin,rightwin,genmatrix,projmatrix,counts) )
names(all.expected) <- rownames(estimates)

# all.resids <- with(long.win, lapply( all.expected, function (x) mapply(function(u,v) (u-v)/sqrt(u),x,counts) ) )
all.resids <- with(long.win, lapply( all.expected, function (x) mapply(function(u,v) qnorm(ppois(as.numeric(v),as.numeric(u))), x, counts, SIMPLIFY=FALSE)) )
for (i in seq_along(all.resids)) for (j in seq_along(all.resids[[i]])) { 
        all.resids[[i]][[j]] <- pmin(100, pmax( -100, all.resids[[i]][[j]] ) )
        dim(all.resids[[i]][[j]]) <- with(long.win,  dim(counts[[1]]) ); 
        dimnames(all.resids[[i]][[j]]) <- with(long.win,  dimnames(counts[[1]]) ) 
    }

resid.counts <- with(long.win, lapply( all.resids[['mle']], countmuts, mutpats=long.mutpats[[2]], leftwin=leftwin) )

resid.counts.mle <- with(long.win, mclapply( long.mutpats, function (mutpats) { lapply( all.resids[['mle']], countmuts, mutpats=mutpats, leftwin=leftwin) }, mc.cores=3 ) )
resid.counts.mean <- with(long.win, mclapply( long.mutpats, function (mutpats) { lapply( all.resids[['mean']], countmuts, mutpats=mutpats, leftwin=leftwin) }, mc.cores=3 ) )

resid.table.mle <- do.call( cbind, lapply( resid.counts.mle[[2]], function (x) { x[1,] } ) )
resid.table.mle <- resid.table.mle[order(rowMeans(resid.table.mle)),] 
resid.table.mean <- do.call( cbind, lapply( resid.counts.mean[[2]], function (x) { x[1,] } ) )
resid.table.mean <- resid.table.mean[order(rowMeans(resid.table.mean)),] 

plot( all.resids[['mle']][[1]], all.resids[['mean']][[1]] )

layout(matrix(1:16,nrow=4))
for (k in 2:17) { hist(mrun$batch[200:1000,k],breaks=40,main=names(mrun$final[k])) }

########## POSTER FIG
pdf( file='poster-results.pdf', width=8, height=3, pointsize=10 )
layout( matrix( 
        c( 1,4,
            2,4,
            3,4 ), nrow=3, byrow=TRUE ), widths=c(2,1) )

# pdf(file='est-mutrates.pdf',width=5,height=3,pointsize=10)
mutlabs <- mutnames(mutpats)
# get reverse-complement by eachother
labperm <- c(1,7,2,5,3,4,6,12,8,11,9,10,13)
par(mar=c(0,3,0,0)+.1)
for (k in labperm[1:12]) {
    first <- (k==labperm[1])
    last <- (k==labperm[7])
    adj <- ( match(k,labperm)%%6 )
    new <- ( adj == 1 )
    if (last) { par(mar=c(3,0,0,0)+.1) }
    if (first) { par(mar=c(0,0,3,0)+.1) }
    hist(mrun$batch[200:1000,1+k], breaks=seq(0,.2,length.out=500), col=adjustcolor(rainbow(20)[k],.5), border=NA, xlim=c(0,.1), xaxt='n', xlab='', main=if (!first) '' else 'posterior density of mutation rates', add=!new, ylim=c(0,60), yaxt='n' )
    if (new) { axis(1,tick=TRUE,labels=last) }
    text(mean(mrun$batch[200:1000,1+k],trim=.1)+.003,10+8*adj,labels=mutlabs[k])
}
hist(mrun$batch[200:1000,1+13], breaks=20, col=grey(.5), xlim=c(0,1), yaxt='n', xlab='', main='')
text(quantile(mrun$batch[200:1000,1+13],.8)+.1,20,labels=mutlabs[13])
# dev.off()

# pdf( file='mle-resids-hist.pdf', width=3, height=3, pointsize=10 )
par(mar=c(4,4,1,1)+.1)
hist( pmin(20,all.resids[['mle']][[1]]), breaks=200, xlab='residual (5,3) tuple count', ylab='abundance', main='' )
#  label tail with these:
biginds <- which( all.resids[['mle']][[1]] > 20, arr.ind=TRUE )
bigresids <- data.frame( mus=rownames(all.resids[['mle']][[1]])[biginds[,1]], rat=colnames(all.resids[['mle']][[1]])[biginds[,2]], resid=all.resids[['mle']][[1]][biginds], stringsAsFactors=FALSE )
legend("topright",legend=apply(bigresids[sample.int(nrow(bigresids),10),1:2],1,paste,collapse='->'),box.lty=0,title='overrepresented',cex=.7)
littleinds <- which( all.resids[['mle']][[1]] < (-10), arr.ind=TRUE )
littleresids <- data.frame( mus=rownames(all.resids[['mle']][[1]])[littleinds[,1]], rat=colnames(all.resids[['mle']][[1]])[littleinds[,2]], resid=all.resids[['mle']][[1]][littleinds], stringsAsFactors=FALSE )
littleresids <- littleresids[order(littleresids$resid),]
legend("topleft",legend=apply(littleresids[1:10,1:2],1,paste,collapse='->'),box.lty=0,title='underrepresented',cex=.7)
# dev.off()

dev.off()
##########

plot( 0, 0, xlim=c(.01,.1), ylim=c(0,100) )
for (k in 2:13) { hist(mrun$batch[200:1000,k],breaks=40,main=names(mrun$final[k]), add=TRUE, col=adjustcolor(rainbow(20)[k],.25), border=NA) }


save(list=ls(),file=with(long.win, paste("resids-",shortwin,"-",leftwin,".RData",sep='')))


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
