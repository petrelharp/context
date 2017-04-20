#!/usr/bin/Rscript

library(contextual)
library(contextutils)
library(simcontext)


## Pattern should be:
# root: X
#  root -> an1: X -> X
#  root -> an2: X -> O
#   an1 -> sp1: X -> O
#   an1 -> sp2: X -> O
#   an2 -> sp3: O -> O
#   an2 -> sp4: O -> X

config <- parse.models( treeify.config( read.config("another-big-tree.json") ) )

clades <- c("an1","an2","sp1","sp2","sp3","sp4")
names(clades) <- clades

genmat.list <- lapply( clades, function (clade) {
        with( list2env(config[[clade]]), 
                makegenmatrix(mutpats,selpats,patlen=1,patterns=bases,fixfn=fixfn,mutrates=mutrates,selcoef=selcoef,selfactors=selfactors,bases=bases,Ne=fixfn.params$Ne)
            )
    } )


pf <- function (ds) { do.call(config$fixfn,c(list(ds),config$fixfn.params)) }

expected.genmats <- cbind( # X->X O->X X->O O->O
                    an1 = c(  0,  0,   0,   0),
                    an2 = c(  0,  0, 10*pf(-.01), 0),
                    sp1 = c(  0,  0, 10*pf(.05),   0),
                    sp2 = c(  0,  0,    pf(1),   0),
                    sp3 = c(  0,  0,   0,   0),
                    sp4 = c(  0, 30*pf(.01),   0,   0)
        )

stopifnot( all(which(expected.genmats>0)==which(sapply(genmat.list,as.vector)>1e-8)) )
stopifnot( all.equal(
        as.vector(expected.genmats),
        as.vector(sapply(genmat.list,as.vector))
    ) )

simseqs <- simseq.tree(100,config)

shortpats <- getpatterns(1,bases=config$bases)

shorts <- lapply( clades, function (clade) {
            counttrans(shortpats, shortpats, simseqs=simseqs[[clade]] )
        } )

expecteds <- list( # X->X O->X X->O O->O
              an1 = c(100,  0,   0,   0),
              an2 = c(  0,  0, 100,   0),
              sp1 = c(  0,  0, 100,   0),
              sp2 = c(  0,  0, 100,   0),
              sp3 = c(  0,  0,   0, 100),
              sp4 = c(  0,100,   0,   0)
        )

stopifnot( all.equal( 
        sapply( lapply( shorts, counts ), as.vector ) ,
        do.call(cbind,expecteds)
    ) )
