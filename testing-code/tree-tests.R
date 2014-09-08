#!/usr/bin/Rscript 

source("../context-inference-fns.R")

config <- treeify.config( read.config("big-tree-model.json") )

genmatrices <- lapply( selfname(c("forward","reverse")), function (x) {
        makegenmatrix( mutpats=config[[x]]$mutpats, bases=config$bases, mutrates=config[[x]]$mutrates, patlen=2, fixfn=null.fixfn )
    } )
genmatrices <- genmatrices[ unlist(config[nodenames(config$tree)]) ]

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=0, fpatterns=getpatterns(1,config$bases) )
root.distrn <- runif( nrow(projmatrix) )
root.distrn <- root.distrn / sum(root.distrn)

transmat <- peel.transmat( config$tree, rowtaxon="sp1", coltaxa=c("sp2","sp3","sp4"), genmatrices=genmatrices, projmatrix=projmatrix, root.distrn=root.distrn, return.list=TRUE )

