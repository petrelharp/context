#!/usr/bin/Rscript 

source("../context-inference-fns.R")

config <- treeify.config( read.config("big-tree-model.json") )

rowtaxa <- "sp1"
coltaxa <- c("sp2","sp3","sp4")

root.node <- get.root(config$tree)
row.node <- match( rowtaxa, nodenames(config$tree) )
col.nodes <- match( coltaxa, nodenames(config$tree) )
# find path from root to rowtaxa:
downpath <- c(row.node)
while( ! root.node %in% downpath ) { downpath <- c( get.parent(downpath[1],config$tree), downpath ) }
# list of active nodes that come off of the path from root to rowtaxa:
up.twigs <- col.nodes[ get.parent(col.nodes,config$tree) %in% downpath ]
# list of nodes in cherries
up.cherries <- setdiff( col.nodes, up.twigs )
while (length(up.cherries)>0) {
}
