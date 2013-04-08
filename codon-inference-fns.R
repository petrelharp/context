#!/usr/bin/R
# NOTE that things like winlen must be defined already, e.g.
source("codons.R")

if (FALSE) {
    # size of window on either side of the focal site
    lwin <- rwin <- 2
    winlen <- lwin+rwin+1

    bases <- c("A","C","G","T")
    nbases <- length(bases)

    patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
    npatt <- nrow(patterns)
    baseind <- nbases^(0:(winlen-1))
    rownames(patterns) <- apply(patterns,1,paste,collapse="")
}

# do stuff mod nbases, essentially.
expandind <- function (i) { i <- i-1; out <- matrix(0,length(i),winlen); for (k in 0:(winlen-1)) { out[,k+1] <- i%%nbases; i <- i%/%nbases }; return(out) }
collapseind <- function (u) { return( u %*% nbases^(0:(winlen-1)) + 1 ) }

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
chind <- function (i,u,k,onebased=1) {
    # index of codon that differs from codon i in substituting u at position k
    if (!is.numeric(u)) { u <- match(u,bases) }
    k <- k-1
    jj <- ( (i-onebased) %/% (nbases^k) ) %% nbases
    ( ( i - onebased - (jj%%nbases)*nbases^k + ((u-1)%%nbases)*nbases^k ) %% nbases^winlen ) + onebased
}
onebase.transitions <- lapply( seq_along(bases), function (u) {
        # list of single-base transitions at the central site
        cbind( 1:npatt, sapply( 1:npatt, chind, u=u, k=lwin+1 ) )
    } )

# row/column indices of one-step transitions
onestep <- lapply( 1:winlen, function (k) {
                header <- if (k==1) { "" } else { apply( do.call( expand.grid, rep( list(bases), k-1L ) ), 1, paste, collapse='' ) }
                footer <- if (k==winlen) { "" } else { apply( do.call( expand.grid, rep( list(bases), winlen-k ) ), 1, paste, collapse='' ) }
                ix <- as.vector( outer( outer( rep(bases,nbases), header, function (x,y) paste(y,x,sep='') ), footer, paste, sep='' ) )
                jx <- as.vector( outer( outer( rep(bases,each=nbases), header, function (x,y) paste(y,x,sep='') ), footer, paste, sep='' ) )
                ii <- as.integer( match( ix, rownames(patterns) ) )
                jj <- as.integer( match( jx, rownames(patterns) ) )
                return( data.frame( i=ii, j=jj ) )
        } )
all.onestep <- do.call( rbind, onestep )

singlebase <- function (baserates) {
    # take vector of off-diagonal entires of rate matrix
    # and return (sparse) matrix of codon transition rates
    # NOTE dgTMatrix indexes ** 0-based **
    rmat <- matrix(0,nbases,nbases)
    rmat[row(rmat)!=col(rmat)] <- baserates
    diag(rmat) <- (-1)*rowSums(rmat)
    with( all.onestep, new( "dgTMatrix", i=i-1L, j=j-1L, x=rep(rmat,winlen*nbases^(winlen-1L)), Dim=c(npatt,npatt) ) )
}

if (FALSE) {
    singlebase.dense <- function (baserates) {
        # take matrix of single-base rates
        # and return matrix of codon transition rates
        # SIMPLE METHOD
        rmat <- matrix(0,nbases,nbases)
        rmat[row(rmat)!=col(rmat)] <- baserates
        diag(rmat) <- (-1)*rowSums(rmat)
        rrmat <- matrix(0, nbases^winlen, nbases^winlen )
        mat.list <- c( list(baserates), rep( list(diag(nbases)), winlen-1L ) )
        for (k in 1:winlen) {
            rrmat <- rrmat + Reduce( kronecker, mat.list )
            mat.list <- c( mat.list[ -1 ], mat.list[1] )
        }
        return(rrmat)
    }

    fullrates <- singlebase(baserates)
    fullrates.dense <- singlebase.dense(baserates)

    # test:
    range( fullrates - singlebase.dense(baserates) )

    # test some
    ijx <- data.frame( ix=sample(rownames(patterns),20), k=sample(1:winlen,20,replace=TRUE), v=sample(bases,20,replace=TRUE), stringsAsFactors=FALSE )
    ijx$jx <- ijx$ix
    ijx$u <- substring(ijx$ix,ijx$k,ijx$k)
    substring(ijx$jx,ijx$k,ijx$k) <- ijx$v
    ijx$i <- match( ijx$ix, rownames(patterns) )
    ijx$j <- match( ijx$jx, rownames(patterns) )
    ijx$rate <- baserates[ as.matrix(ijx[, c("u","v") ]) ]
    ijx$sparse <- fullrates[ as.matrix(ijx[,c("i","j")]) ]
    ijx$dense <- fullrates.dense[ as.matrix(ijx[,c("i","j")]) ]
    ijx$onestep <- sapply( 1:nrow(ijx), function (k) { which( ( all.onestep$i == ijx$i[k] ) & ( all.onestep$j == ijx$j[k] ) ) } )
    ijx <- ijx[ order( abs( ijx$sparse - ijx$rate ) ), ]
}
