# Test inference from simulation
source("codons.R")
source("codon-inference-fns.R")
source("sim-context-fns.R")

# maximum size of pattern
patlen <- 3
lwin <- 0
rwin <- 0

# list of patterns that have mutation rates
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  lapply( combn( bases, 2, simplify=FALSE ), rev),  # single-base rates
        list( c("CG","TG"), c("CG","CA") ), # CpG rates
    NULL )

# mutation rates
mutrates <- runif( length(mutpats) )*1e-8

# list of patterns with selection coefficients
# can be regexps but only using "." and "[...]" (needs regexplen() to work)
selpats <- c(
        "[GC]",
        "[AT]",
        # paste( paste(rep(".",lwin),collapse=''), paste( codons$codon[codons$aa %in% synons], paste(rep(".",rwin),collapse=''), sep='' ), sep='' ),
    NULL )
# selection coefficients
selcoef <- runif( length(selpats) )*1e-4

# check:
stopifnot( patlen >= max( c( sapply(mutpats,function (x) max(nchar(x[1]),nchar(x[2]))), sapply(gsub(".","",selpats,fixed=TRUE),regexplen) ) ) )

# other params
Ne <- 1e4
seqlen <- 100000
tlen <- 6e6

# full instantaneous mutation, and transition matrix
#   RECALL THIS OMITS THE DIAGONAL
genmatrix <- makegenmatrix( mutpats, selpats, patlen=patlen )
genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)

simseqs <- simseq( seqlen, tlen, genmatrix, bases )

# check: number of transitions matches expected?
# normalized prob of transitions
obs.table <- table(simseqs$ntrans[c("i","j")])
obs.p <- with( simseqs, sweep(as.matrix(obs.table),1,rowSums(obs.table),"/") )
exp.p <- sweep(genmatrix,1,max(rowSums(genmatrix)),"/")
diag(exp.p) <- 1-rowSums(exp.p)
pvals <- rep( as.vector(exp.p), as.vector(obs.table) )
stopifnot(all(pvals>0))
stopifnot(all(pvals*length(pvals)>1/100))
# plot(as.vector(exp.p),as.vector(obs.p-exp.p))

## transition probabilities?
# size of window on either side of the focal site
lwin <- 2
rwin <- 2
win <- 1
winlen <- lwin+win+rwin

subtransmatrix <- gettransmatrix(mutpats, mutrates, selpats, selcoef, Ne, tlen, win, lwin, rwin)
counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs, lwin )

# note this is only SORT OF expected
expected <- (seqlen-winlen+1) * (1/nbases)^winlen * subtransmatrix
in.expected <- rowSums(counts) * subtransmatrix
changed <- whichchanged(subtransmatrix,lwin=lwin,win=win)

# looks more or less right
layout(matrix(1:4,2))
plot( as.vector(counts), as.vector(expected), col=1+changed ); abline(0,1)
plot( as.vector(counts), as.vector(in.expected), col=1+changed ); abline(0,1)
plot( as.vector(counts)[changed], as.vector(expected)[changed], col=2 ); abline(0,1)
plot( as.vector(counts)[changed], as.vector(in.expected)[changed], col=2 ); abline(0,1)

# how many events per window?
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=4*winlen))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=2*winlen))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=winlen))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=1))


#######
# testing: number of subsequences

xseqlen <- 1000
xpats <- getpatterns(2)
x <- replicate( 1000, { 
            aseq <- paste( sample(bases,xseqlen,replace=TRUE), collapse="" )
            xmatches <- sapply( lapply( xpats, function (p) gregexpr(paste("(?=",p,")",sep=''),aseq,perl=TRUE) ), "[[", 1 )
            sapply( xmatches, length )
        } )
rownames(x) <- xpats
poisx <- rpois(10000,lambda=xseqlen/length(xpats))
hx <- hist(c(x,poisx), freq=FALSE, plot=FALSE)
histx <- invisible( lapply( seq_along(xpats), function (k) hist( x[k,], col=adjustcolor(rainbow(64)[k],.2), plot=FALSE, breaks=hx$breaks, freq=FALSE ) ) )
matplot( sapply(histx, "[[", "density" ), type='l', lty=c(1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,1), col=1:16 )
legend("topright",legend=xpats,lty=c(1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,1),col=1:16)
lines( hist( poisx, breaks=hx$breaks, plot=FALSE )$density, lwd=2 )
