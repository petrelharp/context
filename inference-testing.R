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

# Compare observed and expected counts:
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

# how many events per window?
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=4*winlen))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=2*winlen))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=winlen))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=1))


#####
# Now, evolve from a common root.

seqlen <- 1e5
tlen <- 6e6
Ne <- 1e4
initseq <- rinitseq(seqlen,bases)
simseqs <- lapply(1:2, function (k) simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef ) )

## transition probabilities?
# size of window on either side of the focal site
lwin <- 2
rwin <- 2
win <- 1
winlen <- lwin+win+rwin

subtransmatrix <- gettransmatrix(mutpats, mutrates, selpats, selcoef, Ne, tlen, win, lwin, rwin)
for (k in seq_along(simseqs)) {
    simseqs[[k]]$counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs[[k]], lwin )
}


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
