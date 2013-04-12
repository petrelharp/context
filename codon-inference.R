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
mutrates <- runif( length(mutmats) )*1e-8
# function to transfer these to list of values in mutation matrix
nmutswitches <- sapply(mutmats,NROW)
whichmut <- function (mutrates) { rep(mutrates,times=nmutswitches) }

# list of patterns with selection coefficients
# can be regexps but only using "." and "[...]" (needs regexplen() to work)
selpats <- c(
        "[GC]",
        "[AT]",
        paste( paste(rep(".",lwin),collapse=''), paste( codons$codon[codons$aa %in% synons], paste(rep(".",rwin),collapse=''), sep='' ), sep='' ),
    NULL )
# number of times each selpat matches each pattern
selmatches <- do.call( rbind, lapply(selpats, function (x) {
        rowSums( sapply( 0:(winlen-regexplen(x)), function (k) {
                    xx <- paste( c(rep(".",k), x, rep(".", winlen-regexplen(x)-k)), collapse='' )
                    grepl( xx, rownames(patterns) ) 
            } ) )
    } ) )

# selection coefficients
selcoef <- runif( length(selpats) )*1e-4
# function to transfer these to list of values (signed) in transition matrix
# these are ( transitions ) x ( selpats ) matrix
fromsel <- selmatches[,  do.call( rbind, mutmats )$i ]
tosel <- selmatches[,  do.call( rbind, mutmats )$j ]
seltrans <- tosel - fromsel
seldiff <- function (selcoef) { as.vector( crossprod(selcoef,seltrans)) }

# Other params
tNe <- (6e6/20)*1e4

# full instantaneous transition matrix
transmatrix <- with( do.call( rbind, mutmats ), new( "dgTMatrix", i=i-1L, j=j-1L, x=whichmut(mutrates)*(1+2*tNe*seldiff(selcoef)), Dim=c(npatt,npatt) ) )

# system.time( replicate( 1, expm( transmatrix, method="Higham08" ) ) )
# with winlen=6 is 114 sec
# with winlen=5 is 2.7 sec
# with winlen=4 is .075 sec
