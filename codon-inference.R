#!/usr/bin/R
source("codons.R")
require(expm)

# size of window on either side of the focal site
lwin <- 0
rwin <- 2
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

mutmatrix <- with( do.call( rbind, mutmats ), new( "dgTMatrix", i=i-1L, j=j-1L, x=rep(mutrates,times=sapply(mutmats,NROW)), Dim=c(npatt,npatt) ) )

system.time( replicate( 50, expm( mutmatrix, method="Higham08" ) ) )

# list of patterns with selection coefficients
selpats <- c(
        codons$codon[codons$aa %in% synons],
    NULL )

# 
