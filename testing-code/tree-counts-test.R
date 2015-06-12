#!/usr/bin/Rscript 

source("../context-inference-fns.R",chdir=TRUE)
source("../sim-context-fns.R",chdir=TRUE)

config <- parse.models( treeify.config( read.config("big-tree-model.json") ) )

simseqs <- simseq.tree(100,config)

longpats <- getpatterns(2,bases=config$bases)
shortpats <- getpatterns(1,bases=config$bases)

##########
# check reordering
longclade <- "sp1"
shortclades <- setdiff(config$tree$tip.label,longclade)
counts.1 <- counttrans.list( list(longpats,shortpats)[c(1,rep.int(2,length(shortclades)))], 
    simseqs=simseqs[c(longclade,shortclades)], 
    leftwin=1, bases=config$bases,
    shift=0 )
counts.2 <- counttrans.list( list(longpats,shortpats)[c(1,rep.int(2,length(shortclades)))], 
    simseqs=simseqs[c(longclade,rev(shortclades))], 
    leftwin=1, bases=config$bases,
    shift=0 )
counts.1.2 <- reorder.counts( counts.1, rev(shortclades) )
counts.2.1 <- reorder.counts( counts.2, shortclades )

stopifnot( all.equal( counts.1.2@counts, counts.2@counts ) )
stopifnot( all.equal( counts.1.2@colpatterns, counts.2@colpatterns ) )
stopifnot( all.equal( counts.2.1@counts, counts.1@counts ) )
stopifnot( all.equal( counts.2.1@colpatterns, counts.1@colpatterns ) )
