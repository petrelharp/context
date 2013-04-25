# Simulate context-dependent mutation rate stuff
source("codons.R")
source("codon-inference-fns.R")

seqlen <- 100000
tlen <- 6e6

#######
# Codons

# maximum size of pattern
patlen <- 3
lwin <- 0
rwin <- 0

# all possible patterns
patterns <- getpatterns(patlen)
npatt <- nrow(patterns)

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

# full instantaneous mutation, and transition matrix
# mutmatrix <- with( do.call( rbind, mutmats ), new( "dgTMatrix", i=i-1L, j=j-1L, x=whichmut(mutrates), Dim=c(npatt,npatt) ) )
genmatrix <- makegenmatrix( mutpats, selpats, patterns )
genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)

# max mutation rate
maxrate <- max( rowSums(genmatrix) )
# mean number of possible changes per window 
meanrate <- maxrate * patlen * tlen
# mean length of a run of connected series of changes (upper bound)
1/(1-ppois(q=0,lambda=meanrate,lower.tail=FALSE))

####
# number and locations of possible changes, ordered by time they occur at
n.events <- rpois(1,lambda=maxrate*tlen*(seqlen-patlen+1))
loc.events <- sample(seqlen-patlen+1,n.events,replace=TRUE)

# count transitions, for debugging
ntrans <- Matrix(0,nrow=nrow(patterns),ncol=nrow(patterns))

# initial and final sequences
finalseq <- initseq <- paste( sample(bases,seqlen,replace=TRUE), collapse="" )
for (k in loc.events) {
    subseq <- substr(finalseq, k, k+patlen-1 )
    msubseq <- match( subseq, rownames(patterns) )
    isubseq <- which( (genmatrix@i + 1) == msubseq )
    # indices of possible replacement patterns
    jsubseq <- (genmatrix@j+1)[isubseq]
    # probabilities of choosing these
    psubseq <- genmatrix@x[isubseq]/maxrate
    replind <- sample( c(msubseq,1+seq_along(patterns)[jsubseq]), 1, prob=c(max(0,1-sum(psubseq)),psubseq) )
    ntrans[msubseq,jsubseq] <- ntrans[msubseq,jsubseq] + 1
    replstr <- c(subseq,rownames(patterns))[replind]
    substr( finalseq, k, k+patlen-1 ) <- replace
}

# check
plot( as.vector( ntrans / sweep(genmatrix,1,maxrate,"/") ) )

## transition probabilities?
# size of window on either side of the focal site
lwin <- 1
rwin <- 1
win <- 3
winlen <- lwin+win+rwin

ipatterns <- getpatterns(winlen)
fpatterns <- getpatterns(win)

# Ok, count occurrences.  Note needs perl "lookahead" to count overlapping patterns.
#   (see http://stackoverflow.com/questions/7878992/finding-the-indexes-of-multiple-overlapping-matching-substrings)
initmatches <- sapply( lapply( rownames(ipatterns), function (p) gregexpr(paste("(?=",p,")",sep=''),initseq,perl=TRUE) ), "[[", 1 )
finalmatches <- sapply( lapply( rownames(fpatterns), function (p) gregexpr(paste("(?=",p,")",sep=''),finalseq,perl=TRUE) ), "[[", 1 )
counts <- Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches) )
for (x in seq_along(initmatches))  { 
    for (y in seq_along(finalmatches)) {
        counts[x,y] <- length(intersect(initmatches[[x]],(-lwin)+finalmatches[[y]]))
    } 
}

stopifnot(sum(counts)==seqlen-winlen+1)

# oops, put this in to makegenmatrix somehow
selpats <- c(
        "[GC]",
        "[AT]",
    NULL )

fullgenmatrix <- makegenmatrix( mutpats, selpats, ipatterns )
fullgenmatrix@x <- update(fullgenmatrix,mutrates,selcoef,Ne)

transmatrix <- expm( tlen * (fullgenmatrix-Diagonal(nrow(fullgenmatrix),rowSums(fullgenmatrix))), method="Higham08" )

subgenmatrix <- collapsepatmatrix( fullgenmatrix, lwin=lwin, rwin=rwin )
subtransmatrix <- transmatrix %*% subgenmatrix
changed <- outer( rownames(ipatterns), rownames(fpatterns), function (x,y) { substr(x,lwin+1,lwin+win)!=y } )

expected <- (seqlen-winlen+1) * (1/nbases)^winlen * subtransmatrix
in.expected <- rowSums(counts) * subtransmatrix

plot( as.vector(counts), as.vector(expected), col=1+changed )
abline(0,1)
plot( as.vector(counts), as.vector(in.expected), col=1+changed )
abline(0,1)



#######
# testing: number of subsequences

xseqlen <- 1000
xpats <- rownames(getpatterns(2))
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
