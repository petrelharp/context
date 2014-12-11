#!/usr/bin/Rscript 

source("../context-inference-fns.R")


###
# Little tree

# config <- treeify.config( read.config("little-tree-model.json") )
json.config <- '
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
    "sp2" : "forward"
}
'
config <- fromJSON(json.config,simplifyMatrix=FALSE)

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )
models <- unlist( config[ nodenames(config$tree) ] )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- runif( nrow(projmatrix) )
root.distrn <- root.distrn / sum(root.distrn)

transmat <- peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )

# check if probabilities sum to 1:
stopifnot( all( ( abs( sapply( transmat, sum ) - 1 ) < sqrt(.Machine$double.eps) ) | ( sapply( transmat, function (x) all( abs(rowSums(x)-1) < sqrt(.Machine$double.eps) ) ) ) ) )


###
# Bigger tree

config <- treeify.config( read.config("big-tree-model.json") )

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )
models <- unlist( config[ nodenames(config$tree) ] )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- runif( nrow(projmatrix) )
root.distrn <- root.distrn / sum(root.distrn)

transmat <- peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3","sp4","sp5","sp6"), models=models, genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )

# check if probabilities sum to 1:
stopifnot( all( ( abs( sapply( transmat, sum ) - 1 ) < sqrt(.Machine$double.eps) ) | ( sapply( transmat, function (x) all( abs(rowSums(x)-1) < sqrt(.Machine$double.eps) ) ) ) ) )


