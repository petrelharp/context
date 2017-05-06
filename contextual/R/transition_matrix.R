
#' Find the Transition Matrix
#'
#' Compute the product of exp(tlen*genmatrix) and projmatrix, either on the
#' left or the right (as `transpose` is TRUE or FALSE) either after a fixed time:
#' exp(tlen*genmatrix) or after a gamma-distributed time
#'
#' @param genmatrix A square sparse matrix.
#' @param projmatrix Another sparse matrix, with number of rows equal to number of rows of genmatrix.
#' @param tlen The amount of time (or the scale, in the gamma distribution).
#' @param shape The shape parameter in the gamma distribution.
#' @param names Attach row and column names to the output?
#' @param transpose Multiply exp(tlen*genmatrix) on the left by t(projmatrix)?
#' @param time Either `fixed` or `gamma`.
#'
#' @export
computetransmatrix <- function( genmatrix, projmatrix, tlen=1, shape=1, names=FALSE, transpose=FALSE, time="fixed",...) {
    if (time=="gamma") {
        if (is.null(dim(projmatrix))) { dim(projmatrix) <- c(length(projmatrix),1) }
        subtransmatrix <- gammaAtv( A=( if (transpose) { Matrix::t(genmatrix) } else {genmatrix} ), scale=tlen, shape=shape, v=projmatrix, ... )
    } else {
        totalrates <- Matrix::rowSums(genmatrix)
        scale.t <- mean(totalrates)
        A <- (1/scale.t) * ( ( if (transpose) { Matrix::t(genmatrix) } else {genmatrix} ) - Diagonal(nrow(genmatrix),totalrates) )
        if (is.null(dim(projmatrix))) { dim(projmatrix) <- c(length(projmatrix),1) }
        subtransmatrix <- sapply( 1:ncol(projmatrix), function (k) { expm::expAtv( A=A, t=tlen*scale.t, v=projmatrix[,k], ... )$eAtv } )
    }
    if (names) {
        rownames(subtransmatrix) <- rownames(genmatrix)
        colnames(subtransmatrix) <- colnames(projmatrix)
    }
    return( subtransmatrix )
}



#' Finds the matrix of conditional probabilities of y given x across a two-branch tree.
#'
#' arguments are lists of two: first the "up" branch (leading from simpler summaries), second the "down"
#'
#' @return Matrix with entry [x,y] the probability of seeing y on the "down
#'   branch" given x was seen on the "up" branch.
getupdowntrans <- function ( genmatrix, projmatrix, mutrates, selcoef,
                            initfreqs, tlens=c(1,1), ... ) {
    otherparams <- list(...)
    op.1 <- lapply( otherparams, "[[", 1 )
    op.2 <- lapply( otherparams, "[[", 2 )
    genmatrix.up <- genmatrix
    genmatrix.up@x <- do.call( update_x, c( list( G=genmatrix, mutrates=mutrates[[1]]*tlens[1], selcoef=selcoef[[1]] ), op.1 ) )
    genmatrix.down <- genmatrix
    genmatrix.down@x <- do.call( update_x, c( list( G=genmatrix, mutrates=mutrates[[2]]*tlens[2], selcoef=selcoef[[2]] ), op.2 ) )
    upbranch <- initfreqs * computetransmatrix( genmatrix.up, projmatrix )   #  prob of root, y
    downbranch <- computetransmatrix( genmatrix.down, initfreqs, transpose=TRUE )  # marginal prob of x
    return( computetransmatrix( genmatrix.down, upbranch, transpose=TRUE ) / as.vector(downbranch) )  # conditional prob of y given x
}

#' "prune" a dangling edge in the "up" direction --
#'    return diag(rootfreqs) * e^( tlen * genmatrix(mut,sel) ) %*% tipmatrix .
#' So, result has
#'   upbranch[i,j] = rootfreqs[i] * E[ tipmatrix[X(tlen),j] | X(0) = i ]
upbranch <- function ( genmatrix, tipmatrix, rootfreqs, mutrates, selcoef, tlen, ... ) {
    if (!missing(mutrates)) { genmatrix@x <- update_x( G=genmatrix, mutrates=mutrates*tlen, selcoef=selcoef, ... ) }
    upbranch <- rootfreqs * computetransmatrix( genmatrix, projmatrix )   #  prob of root, y
    return( upbranch )
}

#' Collapse down branches, rather than up -- like upbranch, but transposed --
#'
#' @return e^( tlen * t(genmatrix(mut,sel)) ) %*% rootmatrix
downbranch <- function ( genmatrix, rootmatrix, mutrates, selcoef, tlen, ... ) {
    if (!missing(mutrates)) { genmatrix@x <- update_x( G=genmatrix, mutrates=mutrates*tlen, selcoef=selcoef, ... ) }
    downbranch <- computetransmatrix( genmatrix.down, rootmatrix, transpose=TRUE )  # marginal prob of x
    return(downbranch)
}


##
# pruning-ish algorithm

cherry.transmats <- function (m1,m2,do.names=FALSE) {
    mm <- m1[,rep(1:ncol(m1),ncol(m2))] * m2[,rep(1:ncol(m2),each=ncol(m1))]
    if (do.names) {
        stopifnot( all( rownames(m1) == rownames(m2) ) )
        rownames(mm) <- rownames(m1)
        colnames(mm) <- paste( colnames(m1)[rep(1:ncol(m1),ncol(m2))], colnames(m2)[rep(1:ncol(m2),each=ncol(m1))], sep=',' )
    }
    return(mm)
}

#' Compute probabilities of all combinations, on a tree
#'
#' Does quite a bit of precomputation to figure out how to (re-)do this computation,
#' but if called with return.list=FALSE will only return the transition matrix
#' that matches with counts (possibly reordered using reorder.counts).
#'
#' @param tree The tree.
#' @param rowtaxon The name of a single tip, or the root - that gets the 'long' patterns
#' @param coltaxa A vector of names of tips - these get 'short' patterns
#' @param modelnames A named character vector of model names, with one name for each node or tip of the tree.
#' @param genmatrices A list of genmatrix'es whose names are model names
#' @param projmatrix A projection matrix that moves from 'long' to 'short' pattern space
#' @param root.distrn A numeric vector of long pattern frequencies at the root
#' @param tlens Lengths of the tree's branches.
#' @param return.list Whether to return everything (default) or just the transition matrices
#' @param debug Whether to keep around some extra information useful for
#'        checking - in particular, row and column names on the resulting matrices.
#'
#' @return A list with the following components:
#'
#' transmats - named by the terminal node of each branch, this gives the transition matrix describing
#'             the probability of seeing the patterns at those tips which this branch separates
#'             from the rowtaxon, given a long pattern at the end of this branch closest to the rowtaxon
#' up  - Gives indices of genmatrices and times for the "up" steps of the peeling algorithm
#' root - Same for the "turning around" step at the root
#' down - Same for moving back down towards the rowtaxon
#' final - Same for the final step connecting to the rowtaxon
#' tlens - Lengths of branches
#' tlens.ord - Unused
#' col.order - 
#' row.node - index of the node that corresponds to rows - transmats[[row.node]]
#'            is the transition matrix we are interested in
#'
#' @details 
#' If `return.list` is TRUE (the default), the output can be used as `setup` to
#' subsequent calls to peel.transmat.compute(), saving on recomputation.
#'
#' @examples
#' # Fast evolution on the branches going to tip3 and the ancestor of (tip1,tip2); but slow on those two branches:
#' tree <- ape::read.tree(text="( (tip1 : 1.0, tip2 : 1.0) node1 : 0.5, tip3 : 1.5 ) root;")
#' G1 <- makegenmatrix( patlen=3, mutpats=list(list(c("XO","OX")),list(c("OX","XO"))), mutrates=c(10,10),
#'                 selpats=list(), selfactors=list(), bases=c("X","O"), fixfn=function (...) { 1 } )
#' G2 <- makegenmatrix( patlen=3, mutpats=list(list(c("XO","OX")),list(c("OX","XO"))), mutrates=c(0.1,0.1),
#'                 selpats=list(), selfactors=list(), bases=c("X","O"), fixfn=function (...) { 1 } )
#' genmatrices <- list(fast=G1, slow=G2)
#' projmatrix <- collapsepatmatrix(rownames(G1), leftwin=1, shortwin=1, bases=c("X","O"))
#' transmat <- peel.transmat(tree, rowtaxon="tip1", coltaxa=c("tip2","tip3"), 
#'                              modelnames=c(tip1="slow", tip2="slow", tip3="fast", node1="fast"), 
#'                              genmatrices=genmatrices, projmatrix=projmatrix,
#'                              root.distrn=rep(1/8,length=8), debug=TRUE, return.list=FALSE )
#' sum(transmat) # equals 1
#'
#' @export
peel.transmat <- function (tree, 
                           rowtaxon, 
                           coltaxa, 
                           modelnames, 
                           genmatrices, 
                           projmatrix, 
                           root.distrn,
                           tlens=tree$edge.length,
                           return.list=TRUE, 
                           debug=FALSE) {
    # indices of nodes:
    setup <- peel.transmat.setup(tree, rowtaxon, coltaxa, modelnames, genmatrices, projmatrix, tlens, debug=debug)
    setup <- peel.transmat.compute( setup, genmatrices, root.distrn, tlens, debug=debug )
    if (return.list) {
        # everything
        return( setup )
    } else {
        # just the answer
        return( setup$transmats[[ setup$row.node ]] )
    }
}

#' Set up for peel.transmat
#'
#' Figures out what order to do things in.
#'
#' @describeIn peel.transmat Set up computation for peel.transmat.compute
#'
#' @export
peel.transmat.setup <- function (tree, 
                                 rowtaxon, 
                                 coltaxa, 
                                 modelnames, 
                                 genmatrices, 
                                 projmatrix, 
                                 tlens=tree$edge.length,
                                 return.list=FALSE, 
                                 debug=FALSE) {
    # indices of nodes:
    root.node <- get.root(tree)
    row.node <- match( rowtaxon, nodenames(tree) )
    col.nodes <- match( coltaxa, nodenames(tree) )
    if (is.na(row.node)) { stop(sprintf("Node %s not found in tree.", rowtaxon)) }
    if (any(is.na(col.nodes))) { stop(sprintf("Nodes %s not found in tree.", paste(coltaxa[is.na(col.nodes)],collapse=" "))) }
    not.root.names <- setdiff(nodenames(tree),rootname(tree))
    if (!all(not.root.names %in% names(modelnames))) { stop("Not all nodes found in names(modelnames).") }
    if (!all(modelnames[not.root.names] %in% names(genmatrices))) { stop("Not all modelnames found in names(genmatrices).") }
    # reorder branch lengths to match terminal node
    tlens.ord <- match(seq_along(nodenames(tree)),tree$edge[,2])
    # long x short transition matrices
    transmats <- lapply( nodenames(tree), function (x) NULL )
    for (x in col.nodes) { transmats[[x]] <- projmatrix }
    # find path from root to rowtaxon:
    downpath <- c(row.node)
    while( ! root.node %in% downpath ) { downpath <- c( get.parent(downpath[1],tree), downpath ) }
    # list of active nodes that come off of the path from root to rowtaxon:
    up.twigs <- col.nodes[ get.parent(col.nodes,tree) %in% downpath ]
    # list of other active nodes
    up.active <- setdiff( col.nodes, up.twigs )
    # .resolve.up, .resolve.down, .resolve.root, and .resolve.final each store seven things:
    #    [1] index of first source in transmats, 
    #    [2] index of second source in transmats, 
    #    [3] index of first genmatrix in genmatrices, 
    #    [4] index of second genmatix in genmatrices, 
    #    [5] index of first tlen in tlens, 
    #    [6] index of second tlen in tlens,
    #    [7] index of destination in transmats
    #
    # This function finds the indices of generator matrices and the indices of tlens
    add.gm <- function (kk) { c( match( modelnames[ nodenames(tree)[kk] ], names(genmatrices) ), tlens.ord[kk] ) }
    # save the order cherries are resolved in
    .resolve.up <- matrix(nrow=0,ncol=7)
    # list of ordering of taxa for columns of each matrix
    col.order <- lapply( nodenames(tree), function (x) c() )
    for (x in col.nodes) { col.order[[x]] <- c(x) }
    ### First, deal with all the upwards branches:
    while (length(up.active)>0) {
        up.cherries <- get.cherries(up.active,tree)
        uppair <- up.cherries[1,]
        # compute and combine transition matrices across each branch
        if (debug) cat( "up: ", uppair, "\n" )
        # transmat1 <- computetransmatrix( genmatrices[[ modelnames[ nodenames(tree)[uppair[1]] ] ]], transmats[[uppair[1]]], tlen=tlens[tlens.ord[uppair[1]]], time="fixed", names=debug)
        # transmat2 <- computetransmatrix( genmatrices[[ modelnames[ nodenames(tree)[uppair[2]] ] ]], transmats[[uppair[2]]], tlen=tlens[tlens.ord[uppair[2]]], time="fixed", names=debug)
        # remove cherry
        up.active <- setdiff( up.active, uppair )
        newly.active <- get.parent( uppair[1], tree )
        # transmats[[newly.active]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
        col.order[[newly.active]] <- c( col.order[[ uppair[1] ]], col.order[[ uppair[2] ]] )
        .resolve.up <- rbind( .resolve.up, c( uppair, add.gm(uppair), newly.active ) )
        if (get.parent(newly.active,tree) %in% downpath) {
            up.twigs <- c( up.twigs, newly.active )
        } else {
            up.active <- c( up.active, newly.active )
        }
    }
    # reorder up.twigs to match downpath
    up.twigs <- up.twigs[ match( get.parent(up.twigs,tree), downpath ) ]
    ## Now move back down from the root (=downpath[1]) to rowtaxon
    # transmats[[ downpath[1] ]] <- sweep( 
    #         computetransmatrix( genmatrices[[ modelnames[nodenames(tree)[ up.twigs[1] ]] ]], transmats[[ up.twigs[1] ]], tlen=tlens[ tlens.ord[up.twigs[1]] ], time="fixed", names=debug ), 
    #     1, root.distrn, "*" )
    col.order[[ downpath[1] ]] <- col.order[[ up.twigs[1] ]]
    # save ordering things are resolved in
    .resolve.root <- c( NA, up.twigs[1], add.gm(c(NA,up.twigs[1])), downpath[1] )
    .resolve.down <- .resolve.final <- matrix(nrow=0,ncol=7)
    # now move back down
    if (length(downpath)>1) {
        for (k in seq_along(up.twigs)[-1]) {
            if (debug) cat("down: ", downpath[k], up.twigs[k], "\n" )
        #     transmat1 <- computetransmatrix( genmatrices[[ modelnames[ nodenames(tree)[downpath[k]] ] ]], transmats[[downpath[k-1]]], tlen=tlens[ tlens.ord[downpath[k]] ], transpose=TRUE, time="fixed", names=debug )
        #     transmat2 <- computetransmatrix( genmatrices[[ modelnames[ nodenames(tree)[up.twigs[k]] ] ]], transmats[[up.twigs[k]]], tlen=tlens[ tlens.ord[up.twigs[k]] ], time="fixed", names=debug )
        #     transmats[[ downpath[k] ]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
            col.order[[ downpath[k] ]] <- c( col.order[[ downpath[k-1] ]], col.order[[ up.twigs[k] ]] )
            .resolve.down <- rbind( .resolve.down, c( downpath[k-1], up.twigs[k], add.gm(c(downpath[k],up.twigs[k])), downpath[k] ) )
        }
       ## And finally, finish off
        if (debug) cat("final: ", row.node, "\n")
        # transmats[[row.node]] <- computetransmatrix( genmatrices[[ modelnames[nodenames(tree)[row.node]] ]],
        #         transmats[[downpath[length(downpath)-1]]], tlen=tlens[ tlens.ord[row.node] ], transpose=TRUE, time="fixed", names=debug )
        col.order[[ row.node ]] <- col.order[[ downpath[length(downpath)-1] ]]
        .resolve.final <- c( downpath[length(downpath)-1], NA, add.gm(c(row.node,NA)), row.node)
        if (debug) cat("order: ", paste( col.order[[ row.node ]], sep=',' ), "\n" )
    }
    return( list( transmats=transmats, 
                    up=.resolve.up, root=.resolve.root, down=.resolve.down, final=.resolve.final, 
                    tlens=tlens, tlens.ord=tlens.ord, col.order=col.order,
                    row.node=row.node  ) 
              )
}


#' Compute probabilities of all combinations, on a tree
#'
#' @describeIn peel.transmat Compute transition matrices
#'
#' @export
peel.transmat.compute <- function (setup, genmatrices, root.distrn, tlens, debug=FALSE, return.list=TRUE ) {
    # do the computation using stuff already set up
    # deal with all the upwards branches
    setup$tlens <- tlens
    for (k.up in seq(1,length.out=nrow(setup$up))) {
        upinds <- setup$up[k.up,]
        transmat1 <- computetransmatrix( genmatrices[[ upinds[3] ]], setup$transmats[[upinds[1]]], 
                                        tlen=setup$tlens[upinds[5]], time="fixed", names=debug)
        transmat2 <- computetransmatrix( genmatrices[[ upinds[4] ]], setup$transmats[[upinds[2]]], 
                                        tlen=setup$tlens[upinds[6]], time="fixed", names=debug)
        setup$transmats[[upinds[7]]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
    }
    # now move back down from the root (=downpath[1]) to rowtaxon
    setup$transmats[[ setup$root[7] ]] <- sweep( 
            computetransmatrix( genmatrices[[ setup$root[4] ]], setup$transmats[[ setup$root[2] ]], 
                               tlen=setup$tlens[ setup$root[6] ], time="fixed", names=debug )
            , 1, root.distrn, "*" )
    # now move back down
    for (k.down in seq(1,length.out=nrow(setup$down))) {
        downinds <- setup$down[k.down,]
        transmat1 <- computetransmatrix( genmatrices[[ downinds[3] ]], setup$transmats[[downinds[1]]], 
                                        tlen=setup$tlens[downinds[5]], transpose=TRUE, time="fixed", names=debug )
        transmat2 <- computetransmatrix( genmatrices[[ downinds[4] ]], setup$transmats[[downinds[2]]], 
                                        tlen=setup$tlens[downinds[6]], time="fixed", names=debug )
        setup$transmats[[ downinds[7] ]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
    }
    # and finish off
    finalinds <- setup$final
    if (length(finalinds)>0) {
        setup$transmats[[ finalinds[7] ]] <- computetransmatrix( genmatrices[[ finalinds[3] ]],
                setup$transmats[[finalinds[1]]], tlen=setup$tlens[finalinds[5]], transpose=TRUE, time="fixed", names=debug )
    }
    if (return.list) {
        return( setup )
    } else {
        return( setup$transmats[[ setup$row.node ]] )
    }
}

