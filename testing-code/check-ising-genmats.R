#!/usr/bin/Rscript

# Run some tests!!
source("../sim-context-fns.R")
source("../context-inference-fns.R")

bases <- c("X","O")

patlen <- 2
mutpats <- list( 
    list( c("O","X") ),
    list( c("X","O") )
    ) 
mutrates <- c(3,5)
selpats <- list(
        c("OX","XO"),
        c("X")
    )
selfactors <- list(
        c(7,1),
        c(1)
    )
selcoef <- c(-2,1)


fixfn <- function (ds,...) { 1/(1+exp(-ds)) }

###### 
# Can make correct generator matrices?

# XX OX XO OO
selpats.nowrap <- c( XX=2, OX=-1, XO=-1, OO=0 )
ds.nowrap <- (-1)*outer( selpats.nowrap, selpats.nowrap, "-" )
selpats.wrap <- c( XX=2, OX=-3, XO=-3, OO=0 )
ds.wrap <- (-1)*outer( selpats.wrap, selpats.wrap, "-" )
selpats.mean <- c( XX=-1+0-1+2, OX=-1-2-1+1, XO=-1-2-1+1, OO=-1+0-1 )
ds.mean <- (-1)*outer( selpats.mean, selpats.mean, "-" )
selpats.factors <- c( XX=2, OX=-13, XO=-1, OO=0 )
ds.factors <- (-1)*outer( selpats.factors, selpats.factors, "-" )
selpats.wrap.factors <- c( XX=2, OX=-15, XO=-15, OO=0 )
ds.wrap.factors <- (-1)*outer( selpats.wrap.factors, selpats.wrap.factors, "-" )
true.mutrates <- rbind( 
                c( 0, 5, 5, 0 ),
                c( 3, 0, 0, 5 ),
                c( 3, 0, 0, 5 ),
                c( 0, 3, 3, 0 )
            )

checkit <- function (x,y) { stopifnot( all.equal( as.vector(x), as.vector(y) ) ) }

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=rep(0,length(selpats)), mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=fixfn ),
        fixfn(0)*true.mutrates
    )

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=fixfn ),
        fixfn(ds.nowrap)*true.mutrates
    )

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, selfactors=selfactors, mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=fixfn ),
        fixfn(ds.factors)*true.mutrates
    )

checkit( 
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=rep(0,length(selpats)), mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=fixfn ),
        fixfn(0)*true.mutrates
    )

checkit( 
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=fixfn ),
        fixfn(ds.wrap)*true.mutrates
    )

checkit( 
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, selfactors=selfactors, mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=fixfn ),
        fixfn(ds.wrap.factors)*true.mutrates
    )

checkit(
        meangenmatrix( leftwin=1, rightwin=1, mutpats=mutpats, selpats=selpats, selcoef=rep(0,length(selpats)), mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=fixfn ),
        fixfn(0)*true.mutrates
    )

checkit(
        meangenmatrix( leftwin=1, rightwin=1, mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=fixfn ),
        fixfn(ds.mean)*true.mutrates
    )

checkit(
        meangenmatrix( leftwin=1, rightwin=1, mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=fixfn ),
        fixfn(ds.mean)*true.mutrates
    )

cat("\n\n  Generator matrices all good!\n")


