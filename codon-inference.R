#!/usr/bin/R
source("codons.R")  # list of codons, etc 

# size of window on either side of the focal site
lwin <- rwin <- 2
winlen <- lwin+rwin+1

bases <- c("A","C","G","T")
nbases <- length(bases)

patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
npatt <- nrow(patterns)
baseind <- nbases^(0:(winlen-1))
rownames(patterns) <- apply(patterns,1,paste,collapse="")

# do stuff mod nbases, essentially.
expandind <- function (i) { out <- rep(0,winlen); for (k in 0:(winlen-1)) { out[k+1] <- i%%nbases; i <- i%/%nbases }; return(out) }
collapseind <- function (u) { return( u %*% nbases^(0:(winlen-1)) ) }

diffind <- function (i) {
    # return indices of patterns differing from i by one site
    j <- i-1
    iind <- matrix(0,length(i),winlen*(nbases-1))
    for (k in 0:(winlen-1)) {
        jj <- j %% nbases
        iind[,k*(nbases-1)+(1:(nbases-1))] <- ( i - 1 - jj*nbases^k + setdiff(0:(nbases-1),jj)*nbases^k ) %% nbases^winlen + 1
        j <- j %/% nbases
    }
    return( iind )
}
chind <- function (i,u,k) {
    # index of codon that differs from codon i in substituting u at position k
    if (!is.integer(u)) { u <- match(u,bases) }
    k <- k-1
    jj <- (i-1 %/% nbases^k ) %% nbases
    ( ( i - 1 - jj*nbases^k + (u-1)*nbases^k ) %% nbases^winlen ) + 1
}

# list of patterns that have mutation rates
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  # single-base rates
        list( c("CG","TG"), c("CG","CA") ), # CpG rates
        # do.call( c, with( subset(codons, aa%in%synons), tapply( levels(codon)[as.numeric(codon)], droplevels(aa), function(x)combn(x, 2, simplify=FALSE) ) ) ),  # synonymous codon rates 
    NULL )


# list of patterns with selection coefficients
selpats <- c(
        codons$codon[codons$aa %in% synons],
    NULL )
