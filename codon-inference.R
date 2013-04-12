#!/usr/bin/R
source("codons.R")
require(expm)

# size of window on either side of the focal site
lwin <- 0
rwin <- 1
win <- 3
winlen <- lwin+win+rwin

patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
npatt <- nrow(patterns)
baseind <- nbases^(0:(winlen-1))
rownames(patterns) <- apply(patterns,1,paste,collapse="")

source("codon-inference-fns.R")


# list of patterns that have mutation rates
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  # single-base rates
        lapply( combn( bases, 2, simplify=FALSE ), rev),  # single-base rates
        list( c("CG","TG"), c("CG","CA") ), # CpG rates
        # do.call( c, with( subset(codons, aa%in%synons), tapply( levels(codon)[as.numeric(codon)], droplevels(aa), function(x)combn(x, 2, simplify=FALSE) ) ) ),  # synonymous codon rates 
    NULL )

# list of matrices with indices of changes corresponding to mutation patterns above
mutmats <- lapply( mutpats, function (x) {
        do.call( rbind, lapply( 1:(winlen-nchar(x[1])+1), function (k) {
                i <- which( substr( rownames(patterns), k, k+nchar(x[1])-1 ) == x[1] )
                replstr <- rownames(patterns)[i]
                substr( replstr, k, k+nchar(x[1])-1 ) <- x[2]
                j <- match( replstr, rownames(patterns) )
                data.frame( i=i, j=j )
            } ) )
    } )

# mutation rates
mutrates <- rep( 1, length(mutmats) )/10
# function to transfer these to list of values in mutation matrix
nmutswitches <- sapply(mutmats,NROW)
whichmut <- function (mutrates) { rep(mutrates,times=nmutswitches) }

# list of patterns with selection coefficients
# can be regexps but only using "." and "[...]"
selpats <- c(
        "[GC]",
        "[AT]",
        paste( paste(rep(".",lwin),collapse=''), paste( codons$codon[codons$aa %in% synons], paste(rep(".",rwin),collapse=''), sep='' ), sep='' ),
    NULL )
regexplen <- function (xx) {
    # length of the matching string... grrr.
    sapply( xx, function (x) {
        y <- diff( c(0,grep("[]\\[]",strsplit(x,"")[[1]],value=FALSE),nchar(x)+1) ) - 1  # lengths of bits in and out of "[]"s
        sum( y[(1 ==  (1:length(y))%%2)] ) + (length(y)-1)/2
    } )
}
# number of times each selpat matches each pattern
selmatches <- lapply(selpats, function (x) {
        rowSums( sapply( 0:(winlen-regexplen(x)), function (k) {
                    xx <- paste( c(rep(".",k), x, rep(".", winlen-regexplen(x)-k)), collapse='' )
                    grepl( xx, rownames(patterns) ) 
            } ) )
    } )

# selection coefficients
selcoef <- runif( length(selpats) )
# function to transfer these to list of values (signed) in transition matrix
fromsel <- sapply( selmatches, 
whichsel <- function (selcoef) {
    selcoef%*%(tosel - fromsel)
}

# full instantaneous transition matrix
transmatrix <- 

mutmatrix <- with( do.call( rbind, mutmats ), new( "dgTMatrix", i=i-1L, j=j-1L, x=whichmut(mutrates)+whichsel(selcoef), Dim=c(npatt,npatt) ) )

# system.time( replicate( 50, expm( mutmatrix, method="Higham08" ) ) )
# .27 sec each for winlen = 4
# 6 sec each for winlen = 5
