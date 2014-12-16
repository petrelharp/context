#!/usr/bin/Rscript 

source("../context-inference-fns.R")


####################################
# Little tree
cat("testing basic operations on a little tree:\n")

config <- treeify.config( read.config("little-tree-model.json") )

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )
models <- unlist( config[ nodenames(config$tree) ] )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- runif( nrow(projmatrix) )
root.distrn <- root.distrn / sum(root.distrn)

peel.config <- peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )
transmat <- peel.config$transmat

# check if probabilities sum to 1:
stopifnot( all( ( abs( sapply( transmat, sum ) - 1 ) < sqrt(.Machine$double.eps) ) | ( sapply( transmat, function (x) all( abs(rowSums(x)-1) < sqrt(.Machine$double.eps) ) ) ) ) )


####################################
# Bigger tree
cat("testing basic operations on a bigger tree:\n")

config <- treeify.config( read.config("big-tree-model.json") )

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )
models <- unlist( config[ nodenames(config$tree) ] )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- runif( nrow(projmatrix) )
root.distrn <- root.distrn / sum(root.distrn)

peel.config <- peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3","sp4","sp5","sp6"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )
transmat <- peel.config$transmat

# check if probabilities sum to 1:
stopifnot( all( ( abs( sapply( transmat, sum ) - 1 ) < sqrt(.Machine$double.eps) ) | ( sapply( transmat, function (x) all( abs(rowSums(x)-1) < sqrt(.Machine$double.eps) ) ) ) ) )


####################################
# check out transition probabilities
cat("Testing pattern probabilities on a little tree.\n")

aa <- 0.5
bb <- 1
cc <- 3
dd <- 1/3
t1 <- 0.1
t2 <- 0.2
json.config <- paste( '
{
    "tree" : [ "(sp1 : ', t1, ', sp2 : ', t2, ') root;" ],
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 0.5, 0.5 ],
    "forward" : {
        "mutpats" : [
            [ [ "X", "O" ] ],
            [ [ "O", "X" ] ]
        ],
        "mutrates" : [ ', aa, ', ', bb, ' ]
    },
    "reverse" : {
        "mutpats" : [
            [ [ "X", "O" ] ],
            [ [ "O", "X" ] ]
        ],
        "mutrates" : [ ', cc, ', ', dd,'  ]
    },
    "sp1" : "forward",
    "sp2" : "reverse"
}
', sep='' )
config <- parse.models( treeify.config( fromJSON(json.config,simplifyMatrix=FALSE) ) )
models <- config.dereference( config, nodenames(config$tree) )

longpats <- getpatterns(2,config$bases)
shortpats <- getpatterns(1,config$bases)
projmatrix <- collapsepatmatrix( ipatterns=longpats, leftwin=1, fpatterns=shortpats )
genmatrices <- lapply( config$.models, function (mm) {
        makegenmatrix( mutpats=config[[mm]]$mutpats, selpats=config[[mm]]$selpats, patterns=longpats,
                               boundary="none", bases=config$bases, fixfn=config[[mm]]$fixfn,
                               mutrates=config[[mm]]$mutrates, selcoef=config[[mm]]$selcoef )
    } )

initfreq.index  <- product.index( longpats=longpats, bases=config$bases ) # which base is at each position in each pattern
initfreqs <- config$initfreqs
root.distrn <- get.root.distrn( initfreqs, initfreq.index )
stopifnot( all.equal(sum(root.distrn),1) )

transmat <- peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa="sp2", models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, debug=TRUE, return.list=FALSE )

stopifnot( all.equal(sum(transmat), 1) )

# truth
fwd.genmat <- matrix( c( # XX  OX   XO  OO
                0,  aa,  aa,   0,  # XX
               bb,   0,   0,  aa,  # OX
               bb,   0,   0,  aa,  # XO
                0,  bb,  bb,   0), # OO
            nrow=4, byrow=TRUE )
rev.genmat <- matrix( c( # XX  OX   XO  OO
                0,  cc,  cc,   0,  # XX
               dd,   0,   0,  cc,  # OX
               dd,   0,   0,  cc,  # XO
                0,  dd,  dd,   0), # OO
            nrow=4, byrow=TRUE )
dimnames(fwd.genmat) <- dimnames(rev.genmat) <- list( longpats, longpats )
require(expm)
fwd.transmat <- expm( 0.1 * (fwd.genmat-diag(rowSums(fwd.genmat))) )
rev.transmat <- expm( 0.2 * (rev.genmat-diag(rowSums(rev.genmat))) )

joint.transmat <- rowSums( sapply( 1:4, function (k) {
            root.distrn[k] * outer( fwd.transmat[k,], rev.transmat[k,] )
        } ) )
dim(joint.transmat) <- c(4,4)
dimnames(joint.transmat) <- list( longpats, longpats )
stopifnot( all.equal(sum(joint.transmat), 1) )

true.transmat <- joint.transmat %*% projmatrix

stopifnot( all.equal( as.vector(true.transmat), as.vector(transmat) ) )
cat("Transition matrix on two-taxon tree all good.\n")

####################################
# bigger tree
cat("Testing pattern probabilities on a bigger tree.\n")

aa <- 0.5
bb <- 1
cc <- 3
dd <- 1/3
t1 <- 0.1
t2 <- 0.2
t3 <- 0.3
t4 <- 0.4
json.config <- paste( '
{
    "tree" : [ "(sp1 : ', t1, ', ( sp2 : ', t2, ', sp3 : ', t3, ') an1 : ', t4, ' ) root;" ],
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 0.5, 0.5 ],
    "forward" : {
        "mutpats" : [
            [ [ "X", "O" ] ],
            [ [ "O", "X" ] ]
        ],
        "mutrates" : [ ', aa, ', ', bb, ' ]
    },
    "reverse" : {
        "mutpats" : [
            [ [ "X", "O" ] ],
            [ [ "O", "X" ] ]
        ],
        "mutrates" : [ ', cc, ', ', dd,'  ]
    },
    "an1" : "forward",
    "sp1" : "forward",
    "sp2" : "reverse",
    "sp3" : "reverse"
}
', sep='' )
config <- parse.models( treeify.config( fromJSON(json.config,simplifyMatrix=FALSE) ) )
models <- config.dereference( config, nodenames(config$tree) )

longpats <- getpatterns(2,config$bases)
shortpats <- getpatterns(1,config$bases)
projmatrix <- collapsepatmatrix( ipatterns=longpats, leftwin=1, fpatterns=shortpats )
big.projmatrix <- collapsepatmatrix( ipatterns=longpats, leftwin=0, fpatterns=longpats )
genmatrices <- lapply( config$.models, function (mm) {
        makegenmatrix( mutpats=config[[mm]]$mutpats, selpats=config[[mm]]$selpats, patterns=longpats,
                               boundary="none", bases=config$bases, fixfn=config[[mm]]$fixfn,
                               mutrates=config[[mm]]$mutrates, selcoef=config[[mm]]$selcoef )
    } )

initfreq.index  <- product.index( longpats=longpats, bases=config$bases ) # which base is at each position in each pattern
initfreqs <- config$initfreqs
root.distrn <- get.root.distrn( initfreqs, initfreq.index )
names(root.distrn) <- longpats
stopifnot( all.equal(sum(root.distrn),1) )

big.transmat <- peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, genmatrices=genmatrices, projmatrix=big.projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, debug=TRUE, return.list=FALSE )

transmat <- peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, debug=TRUE, return.list=FALSE )

stopifnot( all.equal(sum(transmat), 1) )

# truth
fwd.genmat <- matrix( c( # XX  OX   XO  OO
                0,  aa,  aa,   0,  # XX
               bb,   0,   0,  aa,  # OX
               bb,   0,   0,  aa,  # XO
                0,  bb,  bb,   0), # OO
            nrow=4, byrow=TRUE )
rev.genmat <- matrix( c( # XX  OX   XO  OO
                0,  cc,  cc,   0,  # XX
               dd,   0,   0,  cc,  # OX
               dd,   0,   0,  cc,  # XO
                0,  dd,  dd,   0), # OO
            nrow=4, byrow=TRUE )
dimnames(fwd.genmat) <- dimnames(rev.genmat) <- list( longpats, longpats )
require(expm)
sp1.transmat <- expm( t1 * (fwd.genmat-diag(rowSums(fwd.genmat))) )
sp2.transmat <- expm( t2 * (rev.genmat-diag(rowSums(rev.genmat))) )
sp3.transmat <- expm( t3 * (rev.genmat-diag(rowSums(rev.genmat))) )
an1.transmat <- expm( t4 * (fwd.genmat-diag(rowSums(fwd.genmat))) )

root.translist <- lapply( seq_along(root.distrn), function (k) {
            root.distrn[k] * ( 
                an1.transmat[k,1] * outer( sp2.transmat[1,], sp3.transmat[1,], "*" ) +
                an1.transmat[k,2] * outer( sp2.transmat[2,], sp3.transmat[2,], "*" ) +
                an1.transmat[k,3] * outer( sp2.transmat[3,], sp3.transmat[3,], "*" ) +
                an1.transmat[k,4] * outer( sp2.transmat[4,], sp3.transmat[4,], "*" ) )
        } )
names(root.translist) <- names(root.distrn)
# full probabilities
joint.transmat <- t( apply( sp1.transmat, 2, function (xx) {
        xx[1] * root.translist[[1]] + 
        xx[2] * root.translist[[2]] + 
        xx[3] * root.translist[[3]] + 
        xx[4] * root.translist[[4]] 
    } ) )
colnames(joint.transmat) <- outer(longpats,longpats,paste,sep=',')
stopifnot( all.equal(sum(joint.transmat), 1) )

pp <- kronecker(projmatrix,projmatrix)
rownames(pp) <- outer(longpats,longpats,paste,sep=',')
colnames(pp) <- outer(shortpats,shortpats,paste,sep=',')
true.transmat <- joint.transmat %*% pp

stopifnot( all.equal( as.vector(joint.transmat), as.vector(big.transmat) ) )
stopifnot( all.equal( as.vector(true.transmat), as.vector(transmat) ) )
cat("Transition matrix on three-taxon tree good.\n")


###################################
## more peeling: check can create setup then update it
cat("Checking if can update peeling setup.\n")

big.transmat.setup <- peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, genmatrices=genmatrices, projmatrix=big.projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, debug=TRUE, return.list=TRUE )

big.transmat.setup <- peel.transmat.compute( setup=big.transmat.setup, models=models, genmatrices=genmatrices, root.distrn=root.distrn, tlens=config$tree$edge.length )

big.transmat <- big.transmat.setup$transmats[[ big.transmat.setup$row.node ]]

transmat.setup <- peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, debug=TRUE, return.list=TRUE )

transmat.setup <- peel.transmat.compute( setup=transmat.setup, models=models, genmatrices=genmatrices, root.distrn=root.distrn, tlens=config$tree$edge.length )

transmat <- transmat.setup$transmats[[ transmat.setup$row.node ]]


stopifnot( all.equal( as.vector(joint.transmat), as.vector(big.transmat) ) )
stopifnot( all.equal( as.vector(true.transmat), as.vector(transmat) ) )
cat("Can update peeling setup.\n")

