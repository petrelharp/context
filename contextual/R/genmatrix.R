###
# genmatrix code

#' Generator Matrix Class (genmatrix)
#'
#' genmatrix extends the sparse matrix class, carrying along more information.
#'
#' genmatrix gives the instantaneous rate for going from patterns x -> y
#' It is laid out in classical Markov fashion, with the rows indexing the "from" states
#' (headpats) and columns indexing the "to" states (tailpats)
#'
#' @name genmatrix-class
#' @rdname genmatrix-class
#' @importClassesFrom Matrix dgCMatrix 
#' @exportClass genmatrix
setClass("genmatrix", representation(
                 muttrans="Matrix",
                 seltrans="Matrix",
                 bases="character",
                 mutpats="list",
                 selpats="list",
                 selfactors="list",
                 boundary="character",
                 fixfn="function"),
     contains = "dgCMatrix")


#' Make a Generator Matrix
#'
#' Make the generator matrix G on the specified set of patterns,
#' carrying with it the means to quickly update itself.
#'  DON'T do the diagonal (i.e. it must be left empty), so that the updating is easier.
#'
#' Result is a superclass of column-oriented sparse matrices (see
#' ?"dgCMatrix-class" and ?"CsparseMatrix-class"), for which the following
#' vectors implicitly index the vector @x that actually stores the entries:
#'   @i ; (0-based) row numbers of nonzero entries
#'   @p : (0-based) indices of the first nonzero entries in each column (so
#'         diff(x@p) gives the number of nonzero entries by column)
#'
#' @param mutpats List of mutation Tmer motifs.
#' @param selpats List of selection motifs.
#' @param patlen Length of the patterns that index rows and columns of the matrix.
#' @param patterns List of the patterns that index rows and columns of the matrix.
#' @param bases Alphabet used.
#' @param fixfn Fixation function.
#' @param mutrates Vector of mutation rates corresponding to mutpats.
#' @param selcoef Vector of selection coefficients corresponding to selpats.
#' @param selfactors List of numeric vectors of the same structure as selpats of relative selection coefficients.
#' @param boundary What to do with the boundary of the *patterns* ("none" or "wrap").
#'
#' @return A genmatrix object.
#'
#' @details
#'
#' Here is what happens in makegenmatrix:
#'
#' Let G be the generator matrix, i.e. G[i,j] is the instantaneous transition rate from patterns[i] to patterns[j].
#' Construct this as follows:
#'  For each mutation pattern u, let I(u) be the set of index pairs (i,j) such that patterns[i] -> patterns[j] can be obtained by applying u,
#'    ordered in some way, so I(u) = ( (i_1,j_1), (i_2,j_2), ..., (i_{n(u)},j_{n(u)}) ).
#'  Then let I_k be the concatenation of these lists for all mutation patterns in mutpats[[k]].
#'  Also, concatenate the I_k to form I.
#'  Note that the same pair (i,j) can occur in some I_k twice, or in different I_k,
#'    for instance if { C -> T } occurs at rate mutrates[1], and both { CG -> TG, CCG -> CTG } occur at rate mutrates[2],
#'    and patterns[1] == CCG, and patterns[23] == CTG, then the pair (1,23) would occur once in I_1 and twice in I_2.
#'  Now,
#'    G[i,j] = \sum_k mutrates[k] * ( number of times (i,j) appears in I_k ) .
#'  Let (x_1, ..., x_N) be the nonzero entries of G, in some order.
#'    and let f_m = (i_m,j_m) be the row,column index of x_m in G.
#'  We need to compute P so that x = P %*% mutrates, i.e. the N x length(mutpats) matrix
#'    P[m,k] = ( number of times f_m appears in I_k )
#'  To do this, also note that if we define the N x length(I) matrix
#'    Q[m,l] = 1 if the l-th element of I is equal to f_m, and 0 otherwise
#'  and the the length(I) x length(mutpats) matrix
#'    J[l,k] = 1 if the l-th element of I came from I_k, and 0 otherwise,
#'  then
#'    P = Q %*% J .
#'
#' @examples
#' # XO -> OX at rate 3 and OX -> XO at rate 1:
#' G <- makegenmatrix( patlen=3, mutpats=list(list(c("XO","OX")),list(c("OX","XO"))), mutrates=c(3,1),
#'                 selpats=list(), selfactors=list(), bases=c("X","O"), fixfn=function (...) { 1 } )
#' dim(G)
#' G
#' # with selection
#' G <- makegenmatrix( patlen=3, mutpats=list(list(c("XO", "OX")), list(c("OX", "XO"))), mutrates=c(3,1),
#'                 selpats=list(c("XX","OO"),c("X")), selfactors=list(c(2,1),c(1)), selcoef=c(5,1),
#'                 bases=c("X","O"), fixfn=function (...) { 1 } )
#' dim(G)
#' G
#'
#' @export
makegenmatrix <- function (
        mutpats, 
        selpats=list(), 
        patlen,
        patterns=getpatterns(patlen,bases), 
        bases, 
        fixfn, 
        mutrates=rep(1,length(mutpats)), 
        selcoef=rep(1,length(selpats)), 
        selfactors=lapply(selpats,sapply,function(x)1), 
        boundary="none", ...) {
    if (!is.numeric(patlen)|(missing(patlen)&missing(patterns))) { stop("need patlen or patterns") }
    if (missing(patlen)) { patlen <- nchar(as.character(patterns[1])) }
    if ( (length(selpats)>0 && max(sapply(unlist(selpats),nchar))>patlen) | max(sapply(unlist(mutpats),nchar))>patlen ) { stop("some patterns longer than patlen") }
    # mutmats is a list of matrices, with one matrix for each of the mutpatls in mutpatll describing the induced mutation process on the specified patterns (see getmutmats function def'n)
    mutmats <- getmutmats(mutpats,patterns,boundary=boundary)
    # concatenate these matrices into one big matrix, allmutmats;
    #   also compute nmutswitches to keep track of which rows of allmutmats correspond to which of mutpats.
    # the row,column indices of allmutmats are exactly the nonzero entries of the resulting generator matrix,
    #   with the first nmutswitches[1] obtainable by mutpats[[1]], the next nmutswitches[2] obtainable by mutpats[[2]], etcetera.
    #   In other words, nmutswitches[k] gives the number of changes that mutpats[[k]] can induce.
    # note that some rows in allmutmats may be identical:
    #   this means that there is more than one way to apply mutation patterns to change patterns[i] to patterns[j].
    #   We correct for this later.
    allmutmats <- do.call( rbind, mutmats )
    nmutswitches <- sapply(mutmats,NROW)
    # allmutmats is the concatenation of all nonzero entries of the generator matrix,
    #   but they come ordered by which mutpat they correspond to.
    # dgCord is the order in which they will appear in the @x vector of the sparse column-oriented generator matrix,
    #   i.e. ordered by column first, then row. Thus, any identical entries will appear sequentially in dgCord.
    # If there were no duplicated rows in allmutmats, then the nonzero entries of the generator matrix would be
    #    @x == mutrates[ rep(1:length(mutpats),each=nmutswitches) ][ dgCord ]
    dgCord <- order( allmutmats$j, allmutmats$i )
    # Now compute helper matrices:
    #   without selection, each nonzero entry in the generator matrix is a linear combination of mutrates;
    #   muttrans is this linear transformation,
    #   i.e. that @x == muttrans %*% mutrates .
    # This is easily constructed: in the order of allmutmats, this is a nrow(allmutmats) x length(mutpats) matrix,
    #   with 1's along a "fat, irregular, diagonal", the k-th row having nmutswitches[k] 1's, and the rest 0's.
    # We then use dgCord to permute rows.
    # In other words, muttrans[i,j] = 1 if mutpat[[j]] can produce the patterns change of allmutmats[i,].
    # ( Construct via TMatrix, which stores row-and-column indices via @i and @j, then convert to CMatrix. )
    muttrans <- dgTtodgC( new( "dgTMatrix", i=1:sum(nmutswitches)-1L, j=rep(seq_along(mutrates),times=nmutswitches)-1L, x=rep(1,sum(nmutswitches)), Dim=c(sum(nmutswitches),length(mutrates)) ) )
    muttrans <- muttrans[ dgCord , , drop=FALSE ]
    # selection is similar to mutation;
    #   we have to compute differences in numbers of matches of the selection patterns, weighted by selfactors
    if (length(selpats)>0) {
        selmatches <- getselmatches( selpats, patterns, selfactors, boundary=boundary )
        # transfer selection coefficients to selective differences involved in each mutation
        #    these are ( transitions ) x ( selpats ) matrix
        ##  the following has (fromsel-tosel); combined to reduce memory usage
        # Made this sparse in case selpats is large
        ## alternative 1:
        selmatches <- Matrix::t(selmatches)
        fromsel <- selmatches[ allmutmats$i, , drop=FALSE ]
        tosel <- selmatches[  allmutmats$j, , drop=FALSE ]
        seltrans <- ( tosel - fromsel )
        seltrans <- seltrans[dgCord,,drop=FALSE]
        ## alternative 2:
        # seltrans <- ( selmatches[,  allmutmats$j, drop=FALSE ] - selmatches[,  allmutmats$i, drop=FALSE ] )
        # seltrans <- seltrans[,dgCord,drop=FALSE]
        # seltrans <- Matrix(t( seltrans ))
    } else {
        seltrans <- Matrix(numeric(0),nrow=nrow(muttrans),ncol=0)
    }
    # Now deal with the fact there may be duplicated rows in allmutmats:
    #   dups[k] is TRUE if the k-th row of allmutmats[dgCord,] is equal to the (k-1)-th row (see def'n of dgCord),
    #   meaning that they correspond to the same mutpat.
    # We then translate dups to a projection matrix which projects into the unique-mutpat space.
    dups <- c( FALSE, (diff(allmutmats[dgCord,,drop=FALSE][,1])==0) & (diff(allmutmats[dgCord,,drop=FALSE][,2])==0) )
    dupproj <- new("dgTMatrix",i=cumsum(!dups)-1L,j=seq_along(dups)-1L,x=rep(1,length(dups)),Dim=c(nrow(allmutmats)-sum(dups),nrow(allmutmats)))
    muttrans <- dupproj %*% muttrans
    #  and we want to sum mutation rates but not selection probabilities
    #     ... and note that identical rows in allmutmats have identical seltrans entries
    # re-use dupproj: now takes the mean of duplicated entries
    dupproj@x <- as.vector(1/table(dupproj@i))[dupproj@i+1]
    seltrans <- dupproj %*% seltrans
    # Construct the full instantaneous mutation and transition matrix.
    # The `with` function here opens the allmutmats namespace, defining vectors
    #   `i` and `j`, such that an expression like `i=(i-1L)[!dups]` is not
    #   recursive, but rather delivering the new `i` in terms of the allmutmats' `i`.
    genmatrix <- with( allmutmats[dgCord,,drop=FALSE], new( "genmatrix",
            i=(i-1L)[!dups], # Only include a given index if it's not a duplicate.
            p=sapply(0:length(patterns), function (k) sum(j[!dups]<k+1)),
            x=rep(1,sum(!dups)), # Placeholder to get calculated by update_x.
            Dim=c(length(patterns),length(patterns)),
            muttrans=muttrans,
            seltrans=seltrans,
            mutpats=mutpats,
            selpats=selpats,
            selfactors=selfactors,
            bases=bases,
            fixfn=fixfn,
            boundary=boundary
        ) )
    rownames( genmatrix ) <- colnames( genmatrix ) <- patterns
    genmatrix@x <- update_x( genmatrix, mutrates, selcoef, ... )
    return(genmatrix)
}

#' Update the Entries of a Generator Matrix With New Rates.
#'
#' Calculate the new entries of `x` from muttrans et al.
#' use like: genmatrix@x <- update_x(genmatrix,...)
#'
#' @param G A genmatrix object.
#' @param mutrates A new vector of mutation rates.
#' @param selcoef A new vector of selection coefficients (length zero if not needed).
#' @param ... Additional arguments passed to the fixation function, `G@fixfn(s,...)`
#'
#' @return A numeric vector of the same length as G@x.
#'
#' @details
#' The returned value is computed as:
#'   fixfn(seltrans %*% selcoef) * (muttrans %*% mutrates)
#' which works as `seltrans` and `muttrans` have been precomputed appropriately.
#'
#' @export
update_x <- function (G, mutrates, selcoef, ...) {
    fixprob <- if (length(selcoef)>0) { G@fixfn( as.vector(G@seltrans%*%selcoef), ... ) } else { 1 }
    as.vector( G@muttrans %*% mutrates ) * fixprob
}


#' @describeIn makegenmatrix Create a Generator Matrix Averaged Over Possible Boundary States
#' @export
meangenmatrix <- function (leftwin,rightwin,patlen,...) {
    longpatlen <- patlen+leftwin+rightwin
    genmat <- makegenmatrix(...,patlen=longpatlen)  # this is G
    if (longpatlen == patlen) {  # no need to do anything else...
        return(genmat)
    }
    projmat <- collapsepatmatrix(ipatterns=rownames(genmat),leftwin=leftwin,rightwin=rightwin, bases=genmat@bases)  # this is P
    # Divide entries by column sums, then transpose. Makes M, which is now has rows indexed by
    # the short version and columns the usual one.
    meanmat <- Matrix::t( sweep( projmat, 2, Matrix::colSums(projmat), "/" ) )
    pgenmat <- meanmat %*% genmat %*% projmat   # this is H = M G P
    # ii gives the (zero-based) row indices of nonzero entries of pgenmat:
    ii <- pgenmat@i
    # jj gives the (zero-based) column indices of nonzero entries of pgenmat:
    #   see ?"Csparsematrix-class" for explanation of @p.
    jj <- rep(1:ncol(pgenmat),times=diff(pgenmat@p)) - 1L
    nondiag <- ( ii != jj )
    # construct matrix to project from x values in big dgCMatrix to little one:
    #  note that H_ij = M_i. G  P_.j  = sum_kl M_ik G_kl P_lj
    #  ... and we have constructed G and H so we already know which elements are nonzero
    #    and can use this to find the linear transformation
    #    from nonzero elements of G (genmat) to nonzero elements of H (pgenmat)
    ij.H <- 1L + cbind( i=ii[nondiag], j=jj[nondiag] )
    ij.G <- 1L + cbind( i=genmat@i, j=rep(1:ncol(genmat),times=diff(genmat@p))-1L )
    pnonz <- Matrix( 0, nrow=nrow(ij.H), ncol=nrow(ij.G), sparse=TRUE )
    # for-loop to avoid memory hogging
    for (k in 1:nrow(ij.H)) { pnonz[k,] <-  meanmat[ij.H[k,"i"],ij.G[,"i"]] * projmat[ij.G[,"j"],ij.H[k,"j"]] }
    # pnonz <- t( apply( ij.H, 1, function (ij) { meanmat[ij[1],ij.G[,"i"]] * projmat[ij.G[,"j"],ij[2]] } ) )
    pp <- sapply( 0:ncol(pgenmat), function(k) sum(jj[nondiag]<k) )
    meangenmat <- new( "genmatrix",
            i=ii[nondiag],
            p=pp,
            x=pgenmat@x[nondiag],
            Dim=pgenmat@Dim,
            Dimnames=pgenmat@Dimnames,
            muttrans = (pnonz %*% genmat@muttrans),
            seltrans = (pnonz %*% genmat@seltrans),
            bases=genmat@bases,
            mutpats=genmat@mutpats,
            selpats=genmat@selpats,
            selfactors=genmat@selfactors,
            boundary=genmat@boundary,
            fixfn=genmat@fixfn
        )
    args <- list(...)
    if (is.null(args$mutrates)) { args$mutrates <- rep(1,length(genmat@mutpats)) }
    if (is.null(args$selcoef)) { args$selcoef <- rep(1,length(genmat@selpats)) }
    meangenmat@x <- do.call(update_x, c( list(G=meangenmat), args ) )
    return( meangenmat )
}

