#!/usr/bin/R
require(Matrix)
require(expm)
source("expm-simple.R")

getpatterns <- function(winlen) {
    patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
    return( apply(patterns,1,paste,collapse="") )
}

getmutmats <- function(mutpats,patterns) {
    # given pairlist of mutation patterns,
    # return list of matrices with indices of changes corresponding to mutation patterns
    #   i.e. if (i,j) is a row of output[[k]], then patterns[j] can be obtained from patterns[i]
    #   by performing the substitution from mutpats[[i]][1] -> mutpats[[i]][2]
    #   at some location within the string.
    winlen <- nchar(patterns[1])
    lapply( mutpats, function (x) {
        patlen <- nchar(x[1])
        do.call( rbind, lapply( 1:(winlen-patlen+1), function (k) {
                i <- which( substr( patterns, k, k+patlen-1 ) == x[1] ) # matches input
                replstr <- patterns[i]
                substr( replstr, k, k+patlen-1 ) <- x[2] # these, substituted
                j <- match( replstr, patterns )  # indices of mutated strings
                data.frame( i=i, j=j )
            } ) )
    } )
}

fixfn <- function (x,Ne) { 
    # total influx of fixation given selection coefficient (s[to] - s[from]) difference x
    if (length(x)==0) { 1 } else { ifelse( x==0, 1, Ne*expm1(-2*x)/expm1(-2*Ne*x) ) } 
}

# genmatrix extends the sparse matrix class, carrying along more information.
setClass("genmatrix", representation(nmutswitches="integer",seltrans="matrix"), contains = "dgTMatrix")

makegenmatrix <- function (mutpats,selpats,patlen=nchar(patterns[1]),patterns=getpatterns(patlen),mutrates=rep(1,length(mutpats)),selcoef=rep(1,length(selpats)),Ne=1) {
    # make the generator matrix, carrying with it the means to quickly update itself.
    #  DON'T do the diagonal, so that the updating is easier.
    # this gives the instantaneous rate for going from patterns x -> y
    if (!is.numeric(patlen)|(missing(patlen)&missing(patterns))) { stop("need patlen or patterns") }
    if ( (length(selpats)>0 && max(sapply(selpats,regexplen))>patlen) | max(sapply(unlist(mutpats),nchar))>patlen ) { stop("some patterns longer than patlen") }
    # list of matrices with indices of changes corresponding to mutation patterns
    mutmats <- getmutmats(mutpats,patterns)
    allmutmats <- do.call( rbind, mutmats )

    # function to transfer these to list of values in mutation matrix
    mutpatlens <- sapply( sapply(mutpats,"[",1), nchar )
    nmutswitches <- sapply(mutmats,NROW)
    whichmut <- function (mutrates) { rep(mutrates,times=nmutswitches) }

    if (length(selpats)>0) {
        # selmatches[i,j] is number of times selpat[i] matches pattern[j]
        selmatches <- do.call( rbind, lapply(selpats, function (x) {
                rowSums( sapply( 0:(patlen-regexplen(x)), function (k) {
                            xx <- paste( c(rep(".",k), x, rep(".", patlen-regexplen(x)-k)), collapse='' )
                            grepl( xx, patterns ) 
                    } ) )
            } ) )
        # function to transfer selection coefficients to selective differences involved in each mutation
        #    these are ( transitions ) x ( mutpats ) matrix
        #     ... make these sparse?
        fromsel <- selmatches[,  allmutmats$i, drop=FALSE ]
        tosel <- selmatches[,  allmutmats$j, drop=FALSE ]
        seltrans <- tosel - fromsel
        seldiff <- function (selcoef) { as.vector( crossprod(selcoef,seltrans)) }
    } else {
        seltrans <- numeric(0)
        dim(seltrans) <- c(0,0)
        seldiff <- function (...) { 0 }
    }

    # full instantaneous mutation, and transition matrix
    genmatrix <- with( allmutmats, new( "genmatrix", i=i-1L, j=j-1L, 
            x=whichmut(mutrates)*fixfn(seldiff(selcoef),Ne), 
            Dim=c(length(patterns),length(patterns)), 
            nmutswitches=nmutswitches, 
            seltrans=seltrans 
        ) )
    rownames( genmatrix ) <- colnames( genmatrix ) <- patterns
    # diag(genmatrix) <- (-1)*rowSums(genmatrix)  # this makes genmatrix a dgCMatrix

    return(genmatrix)
}

update <- function (G, mutrates, selcoef, Ne) {
    # use like: genmatrix@x <- update(genmatrix,...)
    rep(mutrates,times=G@nmutswitches) * fixfn( as.vector(crossprod(selcoef,G@seltrans)), Ne )
}

collapsepatmatrix <- function (ipatterns, lwin=0, rwin=0, fpatterns=getpatterns(nchar(ipatterns[1])-lwin-rwin) ) {
    # returns a (nbases)^k x (nbases)^k-m matrix projection matrix 
    # map patterns onto the shorter patterns obtained by deleting lwin characters at the start and rwin characters at the end
    patlen <- nchar(ipatterns[1])
    win <- patlen - lwin - rwin
    stopifnot(win>0)
    matchpats <- match( substr(ipatterns,lwin+1L,lwin+win), fpatterns )
    matchmatrix <- new( "dgTMatrix", i=(seq_along(ipatterns)-1L), j=(matchpats-1L), x=rep(1,length(ipatterns)), Dim=c(length(ipatterns),length(fpatterns)) )
    rownames(matchmatrix) <- ipatterns
    colnames(matchmatrix) <- fpatterns
    return( matchmatrix )
}

computetransmatrix <- function( genmatrix, tlen, projmatrix, names=TRUE, ... ) {
    A <- tlen*(genmatrix-Diagonal(nrow(genmatrix),rowSums(genmatrix)))
    subtransmatrix <- sapply( 1:ncol(projmatrix), function (k) { expAtv( A=A, v=projmatrix[,k] )$eAtv } )
    if (names) {
        rownames(subtransmatrix) <- rownames(genmatrix)
        colnames(subtransmatrix) <- colnames(projmatrix)
    }
    subtransmatrix
}

gettransmatrix <- function (mutpats, mutrates, selpats, selcoef, Ne, tlen, win, lwin=0, rwin=0, expm=expm.poisson, ... ) {
    # get reduced transition matrix: given (lwin, win, rwin) context, return probability of pattern in win
    #   note: alternative is expm=expm::expm(x,method="Higham08")
    winlen <- lwin+win+rwin
    fullgenmatrix <- makegenmatrix( mutpats, selpats, patlen=winlen)
    fullgenmatrix@x <- update(fullgenmatrix,mutrates,selcoef,Ne)
    projmatrix <- collapsepatmatrix( ipatterns=rownames(fullgenmatrix), lwin=lwin, rwin=rwin )
    subtransmatrix <- computetransmatrix( genmatrix, tlen, projmatrix, names=TRUE )
    # transmatrix <- expm( tlen * (fullgenmatrix-Diagonal(nrow(fullgenmatrix),rowSums(fullgenmatrix))) )  # exponentiate
    # subtransmatrix <- transmatrix %*% projmatrix        # collapse
    return( subtransmatrix )
}

whichchanged <- function (ipatterns,fpatterns,lwin=0,win=nchar(ipatterns[0])) {
    # return indicator corresponding to entries of output of gettransmatrix that have changed
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    return( outer( ipatterns, fpatterns, function (x,y) { substr(x,lwin+1,lwin+win)!=y } ) )
}

# Misc

regexplen <- function (xx) {
    # length of the string matching a regexp that uses only "." and "[...]" (no other special characters!)
    sapply( xx, function (x) {
        y <- diff( c(0,grep("[]\\[]",strsplit(x,"")[[1]],value=FALSE),nchar(x)+1) ) - 1  # lengths of bits in and out of "[]"s
        sum( y[(1 ==  (1:length(y))%%2)] ) + (length(y)-1)/2
    } )
}

# Unused? 
if (FALSE ){

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



    lwin <- 0
    rwin <- 1
    win <- 3
    winlen <- lwin+rwin+win
    patterns <- getpatterns(winlen)
    npatt <- length(patterns)

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
                    ii <- as.integer( match( ix, patterns ) )
                    jj <- as.integer( match( jx, patterns ) )
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
                splitpats <- do.call(rbind, strsplit(patterns,"") )
                u <- uv[1]; v <- uv[2]
                outer( 1:length(patterns), 1:length(patterns), function (x,y) {
                            rowSums( (splitpats[x,]==u) & (splitpats[y,]==v) )
                        } )
            } )
    dim(ntrans) <- c( length(patterns), length(patterns), length(bases)*length(bases) )
    dimnames(ntrans) <- list( patterns, patterns, outer(bases,bases,paste,sep=".") )
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
    ijx <- data.frame( ix=sample((patterns),20), k=sample(1:winlen,20,replace=TRUE), v=sample(bases,20,replace=TRUE), stringsAsFactors=FALSE )
    ijx$jx <- ijx$ix
    ijx$u <- substring(ijx$ix,ijx$k,ijx$k)
    substring(ijx$jx,ijx$k,ijx$k) <- ijx$v
    ijx$i <- match( ijx$ix, (patterns) )
    ijx$j <- match( ijx$jx, (patterns) )
    ijx$rate <- baserates[ as.matrix(ijx[, c("u","v") ]) ]
    ijx$sparse <- fullrates[ as.matrix(ijx[,c("i","j")]) ]
    ijx$dense <- fullrates.dense[ as.matrix(ijx[,c("i","j")]) ]
    ijx$onestep <- sapply( 1:nrow(ijx), function (k) { which( ( all.onestep$i == ijx$i[k] ) & ( all.onestep$j == ijx$j[k] ) ) } )
    ijx <- ijx[ order( abs( ijx$sparse - ijx$rate ) ), ]
}
