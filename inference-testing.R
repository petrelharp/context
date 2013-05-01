# Test inference from simulation
source("codons.R")
source("codon-inference-fns.R")
source("sim-context-fns.R")

# do diagnostic plots?
diagnostics <- interactive()

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
seqlen <- 1e5
# tlen <- 6e6
tlen <- 6e5

# simulate the sequence
simseqs <- simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef )


# full instantaneous mutation, and transition matrix
#   RECALL THIS OMITS THE DIAGONAL
genmatrix <- makegenmatrix( mutpats, selpats, patlen=patlen )
genmatrix@x <- update(genmatrix,mutrates/(patlen-sapply(sapply(mutpats,"[",1),nchar)+1),selcoef,Ne) # see source of simseq for why rescale
# check: number of (infinitesimal) transitions matches expected?
# normalized prob of transitions
obs.table <- table(simseqs$ntrans[c("i","j")])
obs.p <- with( simseqs, sweep(as.matrix(obs.table),1,rowSums(obs.table),"/") )
exp.p <- sweep(genmatrix,1,max(rowSums(genmatrix)),"/")
diag(exp.p) <- 1-rowSums(exp.p)
pvals <- rep( as.vector(exp.p), as.vector(obs.table) )
stopifnot(all(pvals>0))
stopifnot(all(pvals*length(pvals)>1/100))
if (diagnostics) {
    plot(as.vector(exp.p),as.vector(obs.p)-as.vector(exp.p), xlab="expected infinitesimal transition counts", ylab="residuals")
    abline(0,5); abline(0,-5)
}

## transition probabilities?
# size of window on either side of the focal site
lwin <- 2
rwin <- 2
win <- 1
winlen <- lwin+win+rwin

subtransmatrix <- gettransmatrix(mutpats, mutrates, selpats, selcoef, Ne, tlen, win, lwin, rwin)
counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs, lwin )

# averaged
expected <- (seqlen-winlen+1) * (1/nbases)^winlen * subtransmatrix
# accounting for initial sequence
in.expected <- rowSums(counts) * subtransmatrix
# indicator of counts where patterns have changed
changed <- whichchanged(subtransmatrix,lwin=lwin,win=win)

# huh.
layout(matrix(1:4,2))
plot( as.vector(counts), as.vector(expected), col=1+changed ); abline(0,1)
plot( as.vector(counts), as.vector(in.expected), col=1+changed ); abline(0,1)
plot( as.vector(counts)[changed], as.vector(expected)[changed], col=2 ); abline(0,1)
plot( as.vector(counts)[changed], as.vector(in.expected)[changed], col=2 ); abline(0,1)
# should be mixture of poissons -- compare means in bins
usethese <- ( changed & as.vector(in.expected)>.1 )
exp.bins <- cut( as.vector(in.expected)[usethese], 40 )
exp.bin.vals <- tapply( as.vector(in.expected)[usethese], exp.bins, mean )
points( tapply( as.vector(counts)[usethese], exp.bins, mean ), exp.bin.vals, pch=20 )

if (FALSE) {
    # should be mixture of poissons
    changed.hists <- table( counts=factor(as.vector(counts)[changed],levels=0:max(counts[changed])), expected=round(as.vector(expected)[changed],digits=2) )
    changed.pois <- sapply( 1:ncol(changed.hists), function (k) {
                    x <- as.numeric(rownames(changed.hists))
                    colSums(changed.hists)[k] * dpois( x, lambda=as.numeric(colnames(changed.hists)[k]) )
                } )
    unchanged.hists <- table( counts=factor(as.vector(counts)[!changed],levels=min(counts[!changed]):max(counts[!changed])), expected=round(as.vector(expected)[!changed],digits=2) )
    unchanged.pois <- sapply( 1:ncol(unchanged.hists), function (k) {
                    x <- as.numeric(rownames(unchanged.hists))
                    colSums(unchanged.hists)[k] * dpois( x, lambda=as.numeric(colnames(unchanged.hists)[k]) )
                } )
    matplot(changed.hists,type='l',ylim=range(changed.hists[,-1]),lty=2,col=rainbow(ncol(changed.hists)), lwd=2 )
    matlines(changed.pois,col=rainbow(ncol(changed.hists)), lty=1 )
    matplot(unchanged.hists,type='l',ylim=range(unchanged.hists[,-1]),lty=2,col=rainbow(ncol(unchanged.hists)), lwd=2 )
    matlines(unchanged.pois,col=rainbow(ncol(unchanged.hists)), lty=1 )
}


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
