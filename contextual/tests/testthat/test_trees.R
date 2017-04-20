library(contextual)

library(testthat)
library(jsonlite)

set.seed(23)


####################################
# Little tree
context("testing basic operations on a little tree:\n")

little_model <- '
{
    "comment" : "Simple tree, for testing.",
    "tree" : [ "(sp1 : 0.1, sp2 : 0.2) root;" ],
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 0.5, 0.5 ],
    "forward" : {
        "mutpats" : [
            [ [ "X", "O" ] ],
            [ [ "O", "X" ] ]
        ],
        "mutrates" : [ 1, 0.5 ]
    },
    "reverse" : {
        "mutpats" : [
            [ [ "X", "O" ] ],
            [ [ "O", "X" ] ]
        ],
        "mutrates" : [ 0.5, 1 ]
    },
    "sp1" : "forward",
    "sp2" : "reverse"
}
'
config <- treeify.config( read.config(json=little_model) )

# This model has independent mutations at each site.
# Here are the single-site matrices for each edge:
#  patterns are XX OX XO OO
#  recall that the diagonal is omitted
true_genmats <- list( sp1 = rbind(c(-2,1,1,0),
                                  c(0.5,-1.5,0,1),
                                  c(0.5,0,-1.5,1),
                                  c(0,0.5,0.5,-1)),
                      sp2 = rbind(c(-1,0.5,0.5,0),
                                  c(1,-1.5,0,0.5),
                                  c(1,0,-1.5,0.5),
                                  c(0,1,1,-2)))
zero_diag <- function (x) { diag(x) <- rep(0.0,nrow(x)); x }

true_branchmats <- list( sp1 = expm(0.1 * true_genmats[['sp1']]),
                          sp2 = expm(0.2 * true_genmats[['sp2']]) )

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )

test_that("correct generator matrices", {
              expect_equal( as.numeric(zero_diag(true_genmats[["sp1"]])), as.numeric(genmatrices[["forward"]]) ) 
              expect_equal( as.numeric(zero_diag(true_genmats[["sp2"]])), as.numeric(genmatrices[["reverse"]]) ) 
    } )

models <- unlist( config[ nodenames(config$tree) ] )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- c(4,3,2,1)/10

peel.config <- contextual:::peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2"), models=models, genmatrices=genmatrices, 
                                           projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )
transmat <- peel.config$transmats

the_transmat <- transmat[[peel.config$row.node]]

true_transmats <- list( sp2 = projmatrix,
                        root = ( sweep(true_branchmats$sp2, 1, root.distrn, "*") %*% projmatrix),
                        sp1 = t(true_branchmats$sp1) %*% ( sweep(true_branchmats$sp2, 1, root.distrn, "*") %*% projmatrix) )

test_that("check if probabilities sum to 1:", {
    expect_true( all( ( abs( sapply( transmat, sum ) - 1 ) < sqrt(.Machine$double.eps) ) | ( sapply( transmat, function (x) all( abs(rowSums(x)-1) < sqrt(.Machine$double.eps) ) ) ) ) )
})

test_that("check transition matrices agree", {
              expect_equal(transmat$sp2, true_transmats$sp2)
              expect_equal(as.numeric(transmat$root), as.numeric(true_transmats$root))
              expect_equal(as.numeric(transmat$sp1), as.numeric(true_transmats$sp1))
})


####################################
# Bigger tree
context("testing basic operations on a bigger tree:\n")

config <- treeify.config( read.config("big-tree-model.json") )

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )
models <- unlist( config[ nodenames(config$tree) ] )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- runif( nrow(projmatrix) )
root.distrn <- root.distrn / sum(root.distrn)

peel.config <- contextual:::peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3","sp4","sp5","sp6"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )
transmat <- peel.config$transmats

test_that("check if probabilities sum to 1:", {
    expect_true( all( ( abs( sapply( transmat, sum ) - 1 ) < sqrt(.Machine$double.eps) ) | ( sapply( transmat, function (x) all( abs(rowSums(x)-1) < sqrt(.Machine$double.eps) ) ) ) ) )
})



####################################
# check out transition probabilities
context("Testing pattern probabilities on a little tree.\n")
context("Transition matrix on two-taxon tree all good.\n")

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
config <- contextual:::parse.models( treeify.config( jsonlite::fromJSON(json.config,simplifyMatrix=FALSE) ) )
models <- contextual:::config.dereference( config, nodenames(config$tree) )

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

expect_equal(sum(root.distrn), 1)

transmat <- contextual:::peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa="sp2", models=models, genmatrices=genmatrices, projmatrix=projmatrix, 
                                       root.distrn=root.distrn, tlens=config$tree$edge.length, return.list=FALSE )

expect_equal( sum(transmat), 1)

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

expect_equal(sum(joint.transmat), 1)

true.transmat <- joint.transmat %*% projmatrix

expect_equal( as.vector(true.transmat), as.vector(transmat) )

####################################
# bigger tree
context("Testing pattern probabilities on a bigger tree.\n")

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
config <- contextual:::parse.models( treeify.config( jsonlite::fromJSON(json.config,simplifyMatrix=FALSE) ) )
models <- contextual:::config.dereference( config, nodenames(config$tree) )

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

expect_equal(sum(root.distrn),1)

big.transmat <- contextual:::peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, 
                                           genmatrices=genmatrices, projmatrix=big.projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, return.list=FALSE )

transmat <- contextual:::peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, 
                                       genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, return.list=FALSE )

expect_equal(sum(transmat), 1)

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
expect_equal(sum(joint.transmat), 1)

pp <- kronecker(projmatrix,projmatrix)
rownames(pp) <- outer(longpats,longpats,paste,sep=',')
colnames(pp) <- outer(shortpats,shortpats,paste,sep=',')
true.transmat <- joint.transmat %*% pp

test_that("Transition matrix on three-taxon tree:", {
    expect_equal( as.vector(joint.transmat), as.vector(big.transmat) )
    expect_equal( as.vector(true.transmat), as.vector(transmat) )
})


###################################
## more peeling: check can create setup then update it
context("Checking if can update peeling setup.\n")

big.transmat.setup <- contextual:::peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, genmatrices=genmatrices, 
                                                 projmatrix=big.projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, debug=TRUE, return.list=TRUE )

big.transmat.setup <- peel.transmat.compute( setup=big.transmat.setup, models=models, genmatrices=genmatrices, root.distrn=root.distrn, tlens=config$tree$edge.length )

big.transmat <- big.transmat.setup$transmats[[ big.transmat.setup$row.node ]]

transmat.setup <- contextual:::peel.transmat( tree=config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3"), models=models, genmatrices=genmatrices, 
                                             projmatrix=projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, return.list=TRUE )

transmat.setup <- peel.transmat.compute( setup=transmat.setup, models=models, genmatrices=genmatrices, root.distrn=root.distrn, tlens=config$tree$edge.length )

transmat <- transmat.setup$transmats[[ transmat.setup$row.node ]]


test_that("Can update peeling setup:", {
    expect_equal( as.vector(joint.transmat), as.vector(big.transmat) )
    expect_equal( as.vector(true.transmat), as.vector(transmat) )
})

