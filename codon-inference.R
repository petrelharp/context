#!/usr/bin/R

# size of window on either side of the focal site
lwin <- rwin <- 2
winlen <- lwin+rwin+1

bases <- c("A","C","G","T")
nbases <- length(bases)

patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
npatt <- nrow(patterns)
baseind <- nbases^(0:(winlen-1))

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

# indices, for each patterns, of other pattern differing by at most nchanges
nchanges <- 1
closeby <- data.frame( i=1:npatt )
for 
