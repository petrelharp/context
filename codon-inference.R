#!/usr/bin/R
source("codons.R")
source("codon-inference-fns.R")
require(expm)

# size of window on either side of the focal site
lwin <- 0
rwin <- 1
win <- 3
winlen <- lwin+win+rwin

# all possible patterns
patterns <- getpatterns(winlen)
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
        paste( paste(rep(".",lwin),collapse=''), paste( codons$codon[codons$aa %in% synons], paste(rep(".",rwin),collapse=''), sep='' ), sep='' ),
    NULL )
# selection coefficients
selcoef <- runif( length(selpats) )*1e-4

# Other params
Ne <- 1e4
tlen <- (6e6/20)

# full instantaneous mutation, and transition matrix
# mutmatrix <- with( do.call( rbind, mutmats ), new( "dgTMatrix", i=i-1L, j=j-1L, x=whichmut(mutrates), Dim=c(npatt,npatt) ) )
genmatrix <- makegenmatrix( mutpats, selpats, patterns )
genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)

# transition matrix
transmatrix <- expm( tlen*(genmatrix-Diagonal(nrow(genmatrix),rowSums(genmatrix))), method="Higham08" )

# A <- tlen*genmatrix; system.time( replicate( 1, expm( A, method="Higham08" ) ) )
# with winlen=6 is 151 sec
# with winlen=5 is 3.6 sec
# with winlen=4 is .10 sec

