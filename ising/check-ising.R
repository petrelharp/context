#!/usr/bin/Rscript

# Run some tests!!
source("../sim-context-fns.R")
source("../codon-inference-fns.R")

bases <- c("X","O")

patlen <- 2
mutpats <- list( 
    list( c("O","X") ),
    list( c("X","O") )
    ) 
mutrates <- c(3,5)
selpats <- list(
        c("OO","XX"),
        c("OX","XO"),
        c("O"),
        c("X")
    )
selcoef <- c(0,-2,0,1)

fixfn <- function (ds,...) { ifelse( ds==0, 1, 1/(1+exp(-ds)) ) }

# XX OX XO OO
selpats.nowrap <- c( XX=2, OX=-1, XO=-1, OO=0 )
ds.nowrap <- (-1)*outer( selpats.nowrap, selpats.nowrap, "-" )
selpats.wrap <- c( XX=2, OX=-3, XO=-3, OO=0 )
ds.wrap <- (-1)*outer( selpats.wrap, selpats.wrap, "-" )
selpats.mean <- c( XX=-1+0-1+2, OX=-1-2-1+1, XO=-1-2-1+1, OO=-1+0-1 )
ds.mean <- (-1)*outer( selpats.mean, selpats.mean, "-" )
true.mutrates <- rbind( 
                c( 0, 5, 5, 0 ),
                c( 3, 0, 0, 5 ),
                c( 3, 0, 0, 5 ),
                c( 0, 3, 3, 0 )
            )

checkit <- function (x,y) { stopifnot( all.equal( as.vector(x), as.vector(y) ) ) }

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=rep(0,length(selpats)), mutrates=mutrates, patlen=patlen, boundary="none" ),
        true.mutrates
    )

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none" ),
        fixfn(ds.nowrap)*true.mutrates
    )

checkit( 
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=rep(0,length(selpats)), mutrates=mutrates, patlen=patlen, boundary="wrap" ),
        true.mutrates
    )

checkit( 
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap" ),
        fixfn(ds.wrap)*true.mutrates
    )

checkit(
        meangenmatrix( lwin=1, rwin=1, mutpats=mutpats, selpats=rep(0,length(selpats)), selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none" ),
        true.mutrates
    )

checkit(
        meangenmatrix( lwin=1, rwin=1, mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none" ),
        fixfn(ds.mean)*true.mutrates
    )

checkit(
        meangenmatrix( lwin=1, rwin=1, mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap" ),
        fixfn(ds.mean)*true.mutrates
    )

cat("All good!\n")
