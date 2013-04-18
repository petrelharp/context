#!/usr/bin/R
require(Matrix)

getpatterns <- function(winlen) {
    patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
    npatt <- nrow(patterns)
    baseind <- nbases^(0:(winlen-1))
    rownames(patterns) <- apply(patterns,1,paste,collapse="")
    return( patterns )
}

getmutmats <- function(mutpats,patterns) { 
    winlen <- ncol(patterns)
    lapply( mutpats, function (x) {
        # given pairlist of mutation patterns,
        # return list of matrices with indices of changes corresponding to mutation patterns
        do.call( rbind, lapply( 1:(winlen-nchar(x[1])+1), function (k) {
                i <- which( substr( rownames(patterns), k, k+nchar(x[1])-1 ) == x[1] )
                replstr <- rownames(patterns)[i]
                substr( replstr, k, k+nchar(x[1])-1 ) <- x[2]
                j <- match( replstr, rownames(patterns) )
                data.frame( i=i, j=j )
            } ) )
    } )
}
fixfn <- function (x,Ne) { ifelse( x==0, 1, expm1(2*x)/expm1(2*Ne*x) ) }

setClass("genmatrix", representation(nmutswitches="integer",seltrans="matrix",patnames="character"), contains = "dgTMatrix")


makegenmatrix <- function (mutpats,selpats,patterns,mutrates=rep(1,length(mutpats)),selcoef=rep(1,length(selpats)),Ne=1) {
    # make the generator matrix, carrying with it the means to quickly update itself.
    winlen <- ncol(patterns)
    # list of matrices with indices of changes corresponding to mutation patterns above
    mutmats <- getmutmats(mutpats,patterns)
    allmutmats <- do.call( rbind, mutmats )

    # function to transfer these to list of values in mutation matrix
    nmutswitches <- sapply(mutmats,NROW)
    whichmut <- function (mutrates) { rep(mutrates,times=nmutswitches) }

    # number of times each selpat matches each pattern
    selmatches <- do.call( rbind, lapply(selpats, function (x) {
            rowSums( sapply( 0:(winlen-regexplen(x)), function (k) {
                        xx <- paste( c(rep(".",k), x, rep(".", winlen-regexplen(x)-k)), collapse='' )
                        grepl( xx, rownames(patterns) ) 
                } ) )
        } ) )

    # function to transfer these to list of values (signed) in transition matrix
    # these are ( transitions ) x ( selpats ) matrix
    #  ... make these sparse?
    fromsel <- selmatches[,  allmutmats$i ]
    tosel <- selmatches[,  allmutmats$j ]
    seltrans <- tosel - fromsel
    seldiff <- function (selcoef) { as.vector( crossprod(selcoef,seltrans)) }

    # full instantaneous mutation, and transition matrix
    genmatrix <- with( allmutmats, new( "genmatrix", i=i-1L, j=j-1L, x=whichmut(mutrates)*fixfn(seldiff(selcoef),Ne), Dim=c(npatt,npatt), nmutswitches=nmutswitches, seltrans=seltrans, patnames=rownames(patterns) ) )
    # diag(genmatrix) <- (-1)*rowSums(genmatrix)

    return(genmatrix)
}

collapsepatmatrix <- function (tmat, lwin=0, rwin=0, patnames=tmat@patnames) {
    # reduce an (nbases)^k x (nbases)^k matrix
    #   to a (nbases)^k x (nbases)^k-m matrix
    # by summing over possible combinations on the left and right, with m = lwin+rwin
    winlen <- nchar(patnames[1])
    win <- winlen - lwin - rwin
    stopifnot(win>0)
    fpatnames <- rownames(getpatterns(win))
    matchpats <- match( substr(patnames,lwin,lwin+win-1), fpatnames )
    matchmatrix <- new( "dgTMatrix", i=(seq_along(patnames)-1), j=(matchpats-1), x=rep(1,length(patnames)) )
    return( tmat %*% matchmatrix )
}

update <- function (G, mutrates, selcoef, Ne) {
    # use like: genmatrix@x <- update(genmatrix,...)
    rep(mutrates,times=G@nmutswitches) * fixfn( as.vector(crossprod(selcoef,G@seltrans)), Ne)
}

## Older stuff
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

# Misc

regexplen <- function (xx) {
    # length of the string matching a regexp that uses only "." and "[...]" (no other special characters!)
    sapply( xx, function (x) {
        y <- diff( c(0,grep("[]\\[]",strsplit(x,"")[[1]],value=FALSE),nchar(x)+1) ) - 1  # lengths of bits in and out of "[]"s
        sum( y[(1 ==  (1:length(y))%%2)] ) + (length(y)-1)/2
    } )
}

# Unused?  Takes a long time anyhow.
if (FALSE ){

    lwin <- 0
    rwin <- 1
    win <- 3
    winlen <- lwin+rwin+win
    patterns <- getpatterns(winlen)
    npatt <- nrow(patterns)

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

    ntrans <- apply( expand.grid(bases,bases), 1, function (uv) {
                # list of matrices giving the number of u->v transitions involved for each pattern change
                u <- uv[1]; v <- uv[2]
                outer( 1:nrow(patterns), 1:nrow(patterns), function (x,y) {
                            rowSums( (patterns[x,]==u) & (patterns[y,]==v) )
                        } )
            } )
    dim(ntrans) <- c( nrow(patterns), nrow(patterns), length(bases)*length(bases) )
    dimnames(ntrans) <- list( rownames(patterns), rownames(patterns), outer(bases,bases,paste,sep=".") )
    # total number (nonidentical) transitions
    ttrans <- rowSums(ntrans[,,row( dimnames(ntrans)[[3]] ) != col( dimnames(ntrans)[[3]] )],dims=2)
}

# Testing stuff
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
