# Simulate context-dependent mutation rate stuff
source("codons.R")
source("codon-inference-fns.R")

nbases <- 100000
tlen <- 6e6

# maximum size of pattern
patlen <- 3

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
        codons$codon[codons$aa %in% synons],
    NULL )
# selection coefficients
selcoef <- runif( length(selpats) )*1e-4

# check:
stopifnot( patlen >= max( c( sapply(mutpats,function (x) max(nchar(x[1]),nchar(x[2]))), sapply(gsub(".","",selpats,fixed=TRUE),regexplen) ) ) )

# as mutpats, but right-padded to have a common width
fullmutpats <- lapply( mutpats, function (x)  { sapply( gsub(".","",x,fixed=TRUE), function (y) paste( y, paste(rep(".",patlen-regexplen(y)),collapse=''), sep='' ) ) } )
fullselpats <- sapply( gsub(".","",selpats,fixed=TRUE), function (x)  { paste( x, paste(rep(".",patlen-regexplen(x)),collapse=''), sep='' ) } )

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
n.events <- rpois(1,lambda=maxrate*tlen*(nbases-patlen+1))
loc.events <- sample(nbases-patlen+1,n.events,replace=TRUE)

# initial string -- note right padding to allow rightmost bases to match
finalseq <- initseq <- paste( sample(bases,nbases,replace=TRUE), collapse="" )
for (k in loc.events) {
    subseq <- substr(finalseq, k, k+patlen-1 )
    isubseq <- ( (genmatrix@i + 1) == match( subseq, rownames(patterns) ) )
    # indices of possible replacement patterns
    jsubseq <- (genmatrix@j+1)[isubseq]
    # probabilities of choosing these
    psubseq <- genmatrix@x[isubseq]/maxrate
    replace <- sample( c(subseq,rownames(patterns)[jsubseq]), 1, prob=c(max(0,1-sum(psubseq)),psubseq) )
    substr( finalseq, k, k+patlen-1 ) <- replace
}

# size of window on either side of the focal site
lwin <- 0
rwin <- 1
win <- 3
winlen <- lwin+win+rwin

ipatterns <- getpatterns(winlen)
fpatterns <- getpatterns(win)

# Ok, count occurrences.
initmatches <- sapply( lapply( rownames(ipatterns), gregexpr, initseq ), "[[", 1 )
finalmatches <- sapply( lapply( rownames(fpatterns), gregexpr, finalseq ), "[[", 1 )
counts <- Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches) )
for (x in seq_along(initmatches))  for (y in seq_along(finalmatches)) {
    counts[x,y] <- length(intersect(initmatches[[x]][1],finalmatches[[y]][1]))
}

# oops, put this in to makegenmatrix somehow
selpats <- c(
        "[GC]",
        "[AT]",
        paste( paste(rep(".",lwin),collapse=''), paste( codons$codon[codons$aa %in% synons], paste(rep(".",rwin),collapse=''), sep='' ), sep='' ),
    NULL )

fullgenmatrix <- makegenmatrix( mutpats, selpats, ipatterns )
fullgenmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)

subgenmatrix <- collapsepatmatrix( fullgenmatrix, lwin=lwin, rwin=rwin )
