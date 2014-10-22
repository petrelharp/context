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
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=rep(0,length(selpats)), mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=fixfn ),
        fixfn(0)*true.mutrates
    )

checkit( 
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=fixfn ),
        fixfn(ds.wrap)*true.mutrates
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


######
# Simulates correctly?


#  Long sequences:
tlen <- .1
seqlen <- 1000/tlen
mutrates <- c(1,1)
selcoef <- c(-.5,.5)
simseqs <- simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, bases=bases, fixfn=fixfn )

# counts
left.win <- 2; right.win <- 2; short.win <- 2
long.win <- left.win+short.win+right.win
genmatrix <- makegenmatrix( patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, bases=bases, fixfn=fixfn )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=left.win, rightwin=right.win, bases=genmatrix@bases )
counts <- counttrans( rownames(projmatrix), colnames(projmatrix), simseqs=simseqs, leftwin=left.win, bases=genmatrix@bases )
cwin <- 2
subcounts <- projectcounts( counts, countwin=cwin, lcountwin=0, rcountwin=0 )

all.genmats <- list (
            "simple"=makegenmatrix( patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none", bases=bases, fixfn=fixfn ),
            "wrap"=makegenmatrix( mutpats=mutpats, selpats=selpats, patlen=long.win, mutrates=mutrates*tlen, selcoef=selcoef, boundary="wrap", bases=bases, fixfn=fixfn ),
            "mean"=meangenmatrix( leftwin=1, rightwin=1, patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="wrap", bases=bases, fixfn=fixfn )
        )

all.expected <- lapply( all.genmats, function (genmatrix) {
        expected <- predictcounts( longwin=long.win, shortwin=short.win, leftwin=left.win, initcounts=rowSums(counts), mutrates=tlen*mutrates, selcoef=selcoef, genmatrix=genmatrix, projmatrix=projmatrix )
        subexpected <- projectcounts( counts=expected, lcountwin=0, countwin=cwin, rcountwin=0 )
        return(list( expected=expected, subexpected=subexpected ) )
    } )

all.resids <- lapply( all.expected, function (x) {
            list( full = counts@counts-x[["expected"]]@counts,
                sub = subcounts@counts-x[["subexpected"]]@counts
                ) } )

vars <- sapply( all.resids, function (x) as.vector(x[["sub"]]) )
vars <-  sweep(vars,1,rowMeans(vars),"-")^2 

# different generator matrices give same answers?
stopifnot( all( rowMeans(vars) < 1e-3 ) )

# and these agree with truth pretty well?
stopifnot( all( sapply( all.expected, function (x) { ( all( (0.5 - abs( ppois( as.vector(subcounts), lambda=as.vector(x[["subexpected"]]) ) - 0.5 )) > (1/2)*(1-(1-.001)^(1/length(as.vector(subcounts)))) ) ) } ) ) )

cat("\n\n  Simulations all good!\n")

#####
## Inference works?

genmatrix <- meangenmatrix( leftwin=1, rightwin=1, patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none", bases=bases, fixfn=fixfn )

# params are: mutrates*tlen, selcoef 
likfun <- function (params) {
    # params are: mutrates, selcoef
    # this is collapsed transition matrix
    mutrates <- params[1:nmuts(genmatrix)]
    selcoef <- params[nmuts(genmatrix)+1:nsel(genmatrix)]
    genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef)
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
    # return negative log-likelihood 
    ans <- sum( counts@counts * log(subtransmatrix) )
    if (!is.finite(ans)) { print(paste("Warning: non-finite likelihood with params:",paste(params,collapse=' '))) }
    return(ans)
}

initpar <- c( 2*runif( length(mutpats) ), 2*runif( length(selpats) ) ) # random init
truth <- c( mutrates * tlen, selcoef )  # truth
lbs <- c( rep( 1e-6, nmuts(genmatrix) ), rep( -Inf, nsel(genmatrix) ) )
ubs <- c( rep( 1, nmuts(genmatrix) ), rep( 5, nsel(genmatrix) ) )
parscale <- c( rep( mean(mutrates*tlen), nmuts(genmatrix) ), rep(mean(abs(selcoef)), nsel(genmatrix) ) )
baseval <- likfun(truth)
stopifnot( is.finite(baseval) )

cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3, fnscale=(-1)*abs(baseval), parscale=parscale) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3, fnscale=(-1)*abs(baseval), parscale=parscale) )

estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
colnames(estimates) <- c( paste("muttime",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
estimates$likfun <- apply( estimates, 1, likfun )


# look at observed/expected counts
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            predictcounts( longwin=long.win, shortwin=short.win, leftwin=left.win, initcounts=rowSums(counts), mutrates=x[1:nmuts(genmatrix)], selcoef=x[nmuts(genmatrix)+(1:nsel(genmatrix))], genmatrix=genmatrix, projmatrix=projmatrix )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- 2
subcounts <- projectcounts( countwin=cwin, lcountwin=0, rcountwin=0, counts=counts )
all.subexpected <- lapply( all.expected, function (x) { projectcounts( countwin=cwin, lcountwin=0, rcountwin=0, counts=x ) } )


if (interactive()) {
    layout(matrix(1:4,nrow=2))
    cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
    for (k in 1:4) {
        lord <- order( all.subexpected[["truth"]]@counts[,k] )
        plot( subcounts@counts[lord,k], xaxt='n', xlab='', main=colnames(subcounts)[k] )
        axis(1,at=1:nrow(subcounts),labels=rownames(subcounts)[lord],las=3)
        invisible( lapply(seq_along(all.subexpected),function(j) { lines(all.subexpected[[j]]@counts[lord,k],col=cols[j]) } ) )
        legend("topleft",legend=names(all.subexpected),lty=1,col=cols)
    }
}


stopifnot( all( abs( ppois(as.vector(subcounts),lambda=as.vector(all.subexpected$ans)) - 0.5 ) < .99 ) )
stopifnot( all( abs( ppois(as.vector(subcounts),lambda=as.vector(all.subexpected$cheating)) - 0.5 ) < .99 ) )


cat("\n\n  Inference all good!\n")

