#!/usr/bin/R
require(Matrix,warn.conflicts=FALSE)
require(expm,warn.conflicts=FALSE)
# find what dir we're in
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
PATH <- dirname(frame_files[[length(frame_files)]])
source(paste(PATH,"/expm-simple.R",sep=''))

getpatterns <- function(winlen) {
    patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
    return( apply(patterns,1,paste,collapse="") )
}

npatterns <- function (winlen) {
    return( length(bases)^winlen )
}

getmutmats <- function(mutpats,patterns,boundary=c("wrap","none")) {
    # given mutation patterns,
    #   which can be a list of either pairs or lists of pairs,
    # return list of matrices with (1-based) indices of changes corresponding to mutation patterns
    #   i.e. if (i,j) is a row of output[[k]], then patterns[j] can be obtained from patterns[i]
    #   by performing the substitution from mutpats[[i]][1] -> mutpats[[i]][2]
    #   at some location within the string.
    boundary <- match.arg(boundary)
    winlen <- nchar(patterns[1])
    mutpats <- lapply( mutpats, function (x) { if (is.list(x)) { x } else { list(x) } } )
    lapply( mutpats, function (y) {  # y is list of short pattern pairs
            do.call( rbind, lapply(y, function (x) {  # x is short pattern pair (from, to)
                patlen <- nchar(x[1])
                switch( boundary, 
                    wrap={ # patterns are circular
                        wpatterns <- paste( patterns, substr(patterns,1,patlen), sep='' )
                        maxshift <- winlen
                        subsfun <- function (pat,topat,k) { wrapsubstr( pat, k, k+patlen-1 ) <- topat; return(pat) }
                    },
                    none={
                        wpatterns <- patterns
                        maxshift <- winlen-patlen+1
                        subsfun <- function (pat,topat,k) { substr( pat, k, k+patlen-1 ) <- topat; return(pat) }
                    }
                )
                do.call( rbind, lapply( 1:maxshift, function (k) {  # k is position of short pattern in long pattern
                        i <- which( substr( wpatterns, k, k+patlen-1 ) == x[1] ) # which patterns match short from-pattern?
                        replstr <- subsfun( patterns[i], x[2], k ) # substitute in to-pattern
                        j <- match( replstr, patterns )  # indices of mutated strings
                        data.frame( i=i, j=j )
                    } ) )
        } ) )
    } )
}

getselmatches <- function (selpats, patterns, boundary=c("wrap","none"), names=FALSE) {
    # selpats can be a vector or a list of vectors,
    #  each element gets one parameter
    # selmatches[i,j] is number of times anything in selpat[[i]] matches pattern[j]
    boundary <- match.arg(boundary)
    substrfun <- switch( boundary, wrap=wrapsubstr, none=substr )
    patlen <- nchar(patterns[1])
    if (!is.list(selpats)) { selpats <- as.list(selpats) }
    selmatches <- do.call( rbind, lapply(selpats, function (y) {
            rowSums( sapply(y, function (x) {
                    maxshift <- patlen - switch( boundary, wrap=1, none=nchar(x) )
                    rowSums( sapply( 0:maxshift, function (k) {
                                x == substrfun( patterns, 1+k, k+nchar(x) )
                                # xx <- paste( c(rep(".",k), x, rep(".", patlen-regexplen(x)-k)), collapse='' )
                                # grepl( xx, patterns ) 
                        } ) )
                } ) )
        } ) )
    if (names) {
        rownames(selmatches) <- names(selpats)
        colnames(selmatches) <- patterns
    }
    return(selmatches)
}


popgen.fixfn <- function (ds,Ne) { 
    # total influx of fixation given selection coefficient (s[to] - s[from]) difference ds
    if (length(ds)==0) { 1 } else { ifelse( ds==0, 1, Ne*expm1(-2*ds)/expm1(-2*Ne*ds) ) } 
}

# genmatrix extends the sparse matrix class, carrying along more information.
setClass("genmatrix", representation(muttrans="Matrix",seltrans="Matrix",mutrates="numeric",selcoef="numeric"), contains = "dgCMatrix")

makegenmatrix <- function (mutpats,selpats,patlen=nchar(patterns[1]),patterns=getpatterns(patlen),mutrates=rep(1,length(mutpats)),selcoef=rep(1,length(selpats)), boundary=c("wrap","none"), ...) {
    # make the generator matrix, carrying with it the means to quickly update itself.
    #  DON'T do the diagonal, so that the updating is easier.
    # this gives the instantaneous rate for going from patterns x -> y
    if (!is.numeric(patlen)|(missing(patlen)&missing(patterns))) { stop("need patlen or patterns") }
    if ( (length(selpats)>0 && max(sapply(unlist(selpats),nchar))>patlen) | max(sapply(unlist(mutpats),nchar))>patlen ) { stop("some patterns longer than patlen") }
    # list of matrices with indices of changes corresponding to mutation patterns
    mutmats <- getmutmats(mutpats,patterns,boundary=boundary)
    allmutmats <- do.call( rbind, mutmats )
    # convert to dgCMatrix format
    dgCord <- order( allmutmats$j, allmutmats$i )
    # function to transfer these to list of values in mutation matrix
    nmutswitches <- sapply(mutmats,NROW)
    muttrans <- dgTtodgC( new( "dgTMatrix", i=1:sum(nmutswitches)-1L, j=rep(seq_along(mutrates),times=nmutswitches)-1L, x=rep(1,sum(nmutswitches)), Dim=c(sum(nmutswitches),length(mutrates)) ) )
    muttrans <- muttrans[ dgCord , , drop=FALSE ]
    # selection?
    if (length(selpats)>0) {
        selmatches <- getselmatches( selpats, patterns, boundary=boundary )
        # transfer selection coefficients to selective differences involved in each mutation
        #    these are ( transitions ) x ( mutpats ) matrix
        #     ... make these sparse?
        fromsel <- selmatches[,  allmutmats$i, drop=FALSE ]
        tosel <- selmatches[,  allmutmats$j, drop=FALSE ]
        seltrans <- Matrix(t( (tosel - fromsel)[,dgCord,drop=FALSE] ))
    } else {
        seltrans <- Matrix(numeric(0),nrow=nrow(muttrans),ncol=0)
    }
    # there may be duplicated rows (matching multiple patterns); deal with this
    dups <- c(FALSE,diff(allmutmats[dgCord,][,1])==0)
    dupproj <- new("dgTMatrix",i=cumsum(!dups)-1L,j=seq_along(dups)-1L,x=rep(1,length(dups)),Dim=c(nrow(allmutmats)-sum(dups),nrow(allmutmats)))
    muttrans <- dupproj %*% muttrans
    seltrans <- dupproj %*% seltrans
    # full instantaneous mutation, and transition matrix
    genmatrix <- with( allmutmats[dgCord,], new( "genmatrix", 
            i=(i-1L)[!dups], 
            p=sapply(0:length(patterns), function (k) sum(j[!dups]<k+1)),
            x=rep(1,sum(!dups)),
            Dim=c(length(patterns),length(patterns)), 
            muttrans=muttrans,
            seltrans=seltrans,
            mutrates=mutrates,
            selcoef=selcoef
        ) )
    rownames( genmatrix ) <- colnames( genmatrix ) <- patterns
    genmatrix@x <- update( genmatrix, mutrates, selcoef, ... )
    # diag(genmatrix) <- (-1)*rowSums(genmatrix)  # this makes genmatrix a dgCMatrix
    return(genmatrix)
}

update <- function (G, mutrates=G@mutrates, selcoef=G@selcoef, ...) {
    # use like: genmatrix@x <- update(genmatrix,...)
    fixprob <- if (length(selcoef)>0) { fixfn( as.vector(G@seltrans%*%selcoef), ... ) } else { 1 }
    as.vector( G@muttrans %*% mutrates ) * fixprob 
}

collapsepatmatrix <- function (ipatterns, lwin=0, rwin=nchar(ipatterns[1])-nchar(fpatterns[1])-lwin, fpatterns=getpatterns(nchar(ipatterns[1])-lwin-rwin) ) {
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

meangenmatrix <- function (lwin,rwin,patlen,...) {
    # create a generator matrix that averages over possible adjacent states
    longpatlen <- patlen+lwin+rwin
    genmat <- makegenmatrix(...,patlen=longpatlen)  # this is G
    projmat <- collapsepatmatrix(ipatterns=rownames(genmat),lwin=lwin,rwin=rwin)  # this is P
    meanmat <- t( sweep( projmat, 2, colSums(projmat), "/" ) )  # this is M
    pgenmat <- meanmat %*% genmat %*% projmat   # this is H = M G P
    ii <- pgenmat@i
    jj <- rep(1:ncol(pgenmat),times=diff(pgenmat@p)) - 1L
    nondiag <- ( ii != jj )
    # construct matrix to project from x values in big dgCMatrix to little one:
    #  note that H_ij = M_i. G  P_.j  = sum_kl M_ik G_kl P_lj
    #  ... and we have constructed G and H so we know which the nonzero elements are already
    #    and can use this to find the linear transformation 
    #    from nonzero elemetns of G to nonzero elements of H
    ij.H <- 1L + cbind( i=ii[nondiag], j=jj[nondiag] )
    ij.G <- 1L + cbind( i=genmat@i, j=rep(1:ncol(genmat),times=diff(genmat@p))-1L )
    # for-loop to avoid memory hogging
    pnonz <- Matrix( 0, nrow=nrow(ij.H), ncol=nrow(ij.G), sparse=TRUE )
    for (k in 1:nrow(ij.H)) { pnonz[k,] <-  meanmat[ij.H[k,"i"],ij.G[,"i"]] * projmat[ij.G[,"j"],ij.H[k,"j"]] }
    # pnonz <- t( apply( ij.H, 1, function (ij) { meanmat[ij[1],ij.G[,"i"]] * projmat[ij.G[,"j"],ij[2]] } ) )
    pp <- sapply( 0:ncol(pgenmat), function(k) sum(jj[nondiag]<k) )
    meangenmat <- new( "genmatrix", i=ii[nondiag], p=pp, x=pgenmat@x[nondiag], Dim=pgenmat@Dim, Dimnames=pgenmat@Dimnames,
        muttrans = (pnonz %*% genmat@muttrans), seltrans = (pnonz %*% genmat@seltrans),
        mutrates=genmat@mutrates, selcoef=genmat@selcoef )
    meangenmat@x <- update(meangenmat)
    return( meangenmat )
}

computetransmatrix <- function( genmatrix, projmatrix, tlen=1, names=FALSE, transpose=FALSE, ... ) {
    # Compute the product of exp(tlen*genmatrix) and projmatrix, either on the left or the right (as transpose is true or false)
    # tlen confounded with mutation parameters... best to leave that out of here...
    A <- if (tlen==1) { genmatrix - Diagonal(nrow(genmatrix),rowSums(genmatrix)) } else {  tlen*(genmatrix-Diagonal(nrow(genmatrix),rowSums(genmatrix))) }
    if (transpose) { A <- t(A) }
    if (is.null(dim(projmatrix))) { dim(projmatrix) <- c(length(projmatrix),1) }
    subtransmatrix <- sapply( 1:ncol(projmatrix), function (k) { expAtv( A=A, v=projmatrix[,k] )$eAtv } )
    if (names) {
        rownames(subtransmatrix) <- rownames(genmatrix)
        colnames(subtransmatrix) <- colnames(projmatrix)
    }
    subtransmatrix
}

getupdowntrans <- function ( genmatrix, projmatrix, mutrates, selcoef, Ne, initfreqs, tlens=c(1,1) ) {
    # arguments are lists of two: first the "up" branch (leading from simpler summaries), second the "down"
    # returns matrix with entry [x,y] the probability of seeing y on the "down branch" given x was seen on the "up" branch.
    genmatrix.up <- genmatrix
    genmatrix.up@x <- update(genmatrix,mutrates[[1]]*tlens[1],selcoef[[1]],Ne[1])
    genmatrix.down <- genmatrix
    genmatrix.down@x <- update(genmatrix,mutrates[[2]]*tlens[2],selcoef[[2]],Ne[2])
    upbranch <- initfreqs * computetransmatrix( genmatrix.up, projmatrix )   #  prob of root, y
    downbranch <- computetransmatrix( genmatrix.down, initfreqs, transpose=TRUE )  # marginal prob of x
    return( computetransmatrix( genmatrix.down, upbranch, transpose=TRUE ) / as.vector(downbranch) )  # conditional prob of y given x
}


predictcounts <- function (win, lwin=0, rwin=0, initcounts, mutrates, selcoef, mutpats, selpats, genmatrix, projmatrix, ... ) {
    # Compute expected counts of paired patterns:
    winlen <- lwin+win+rwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=lwin+win+rwin,mutpats=mutpats,selpats=selpats, ...) }
    if (!missing(mutrates)) { genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef,...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin ) }
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )
    if (missing(initcounts)) { initcounts <- 1 }
    fullcounts <- initcounts * subtransmatrix
    return( fullcounts )
}

projectcounts <- function( lwin, lcountwin, countwin, rcountwin, counts ) {
    # collapse: valid shift ranges is
    #  (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
    #  using counts for (lwin,win,rwin) compute counts for (lcountwin,countwin,rcountwin).
    if ( max(0L,lwin-lcountwin) > (lwin+win)-(lcountwin+countwin)+min(0L,rwin-rcountwin) ) { stop("unreconcilable windows specified.") }
    winlen <- nchar(rownames(counts)[1])
    win <- nchar(colnames(counts)[1])
    rwin <- winlen-win-lwin
    pcounts <- matrix(0,nrow=npatterns(lcountwin+countwin+rcountwin),ncol=npatterns(countwin))
    for (k in max(0L,lwin-lcountwin):((lwin+win)-(lcountwin+countwin)+min(0L,rwin-rcountwin))) {
        lpmat <- collapsepatmatrix( ipatterns=rownames(counts), lwin=k, rwin=winlen-(k+lcountwin+countwin+rcountwin) )
        rpmat <- collapsepatmatrix( ipatterns=colnames(counts), lwin=k+lcountwin-lwin, rwin=win-(k+lcountwin-lwin+countwin) )
        pcounts <- pcounts + t(lpmat) %*% counts %*% (rpmat)
    }
    dimnames(pcounts) <- list( colnames(lpmat), rownames(rpmat) )
    return(pcounts)
}

whichchanged <- function (ipatterns,fpatterns,lwin=0,win=nchar(ipatterns[0])) {
    # return indicator corresponding to entries of output of gettransmatrix that have changed
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    return( outer( ipatterns, fpatterns, function (x,y) { substr(x,lwin+1,lwin+win)!=y } ) )
}

leftchanged <- function (ipatterns,fpatterns,lwin=0,win=nchar(ipatterns[0])) {
    # return indicator corresponding to whether leftmost element changed
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    return( outer( ipatterns, fpatterns, function (x,y) { substr(x,lwin+1,lwin+1)!=substr(y,1,1) } ) )
}

getlikfun <- function (nmuts,nsel,genmatrix,projmatrix,const=0) {
    return( function (params) {
        # params are: mutrates, selcoef, Ne 
        mutrates <- params[1:nmuts]
        selcoef <- params[nmuts+(1:nsel)]
        Ne <- params[nmuts+nsel+1]
        # tlen <- params[nmuts+nsel+2]  # confounded.
        # this is collapsed transition matrix
        genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
        # return negative log-likelihood 
        (-1) * sum( counts * log(subtransmatrix) ) + const
    } )
}


# Misc

wrapsubstr <- function (x,start,stop) {
    if (all(nchar(x)==0)) { return("") }
    while( any(nchar(x)<stop) ) {
        x <- paste(x,x,sep='')
    }
    substr(x,start,stop)
}

"wrapsubstr<-" <- function (x,start,stop,value) {
    if (length(value)<length(x)) { value <- rep(value,length(x)) }
    xlen <- nchar(x)
    k <- rep(1,length(x))
    while( TRUE ) {
        stop <- ifelse( xlen>0, stop - xlen*((start-1)%/%xlen), 0 )
        start <- ifelse( xlen>0, (start-1)%%xlen+1, 0 )
        thisstop <- pmin(xlen,stop)
        substr(x,start,thisstop) <- substr(value,k,k+thisstop-start)
        k <- k+thisstop-start+1
        start <- thisstop+1
        if( all(stop<start) ) { break; }
    }
    return(x)
}

dgTtodgC <- function (M) {
    # convert between matrix classes.  For understanding.
    ijx <- data.frame( i=M@i, j=M@j, x=M@x )
    ijx <- ijx[ order( ijx$j, ijx$i ), ]
    with(ijx, new( "dgCMatrix", i=i, p=sapply(0:ncol(M), function(k) sum(j<k)), x=x, Dim=dim(M) ) )
}

# Unused? 
if (FALSE ){

regexplen <- function (xx) {
    # length of the string matching a regexp that uses only "." and "[...]" (no other special characters!)
    sapply( xx, function (x) {
        y <- diff( c(0,grep("[]\\[]",strsplit(x,"")[[1]],value=FALSE),nchar(x)+1) ) - 1  # lengths of bits in and out of "[]"s
        sum( y[(1 ==  (1:length(y))%%2)] ) + (length(y)-1)/2
    } )
}

gettransmatrix <- function (mutpats, mutrates, selpats, selcoef, Ne, tlen=1, win, lwin=0, rwin=0, expm=expm.poisson, ... ) {
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
