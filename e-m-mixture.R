#!/usr/bin/R
library(expm)

source("codons.R")  # list of codons, etc 
source("codon-inference-fns.R")

# size of window on either side of the focal site
lwin <- rwin <- 1
winlen <- lwin+rwin+1

bases <- c("A","C","G","T")
nbases <- length(bases)

patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
npatt <- nrow(patterns)
baseind <- nbases^(0:(winlen-1))
rownames(patterns) <- apply(patterns,1,paste,collapse="")

baserates <- runif(nbases*(nbases-1))

# transition matrix for indep't sites
t <- 2
fullrates <- singlebase( baserates )
pmat <- expm( t*fullrates, method="Higham08" )

# prob of selection
p <- .001

# generate data
nsamples <- 1000000
xy <- data.frame( x=sample(1:nrow(patterns),nsamples,replace=TRUE), sel=rbinom(nsamples,size=1,prob=p) )
xy$y <- NA
for (x in unique(xy$x)) {
    xy$y[ xy$x==x ] <- sample( 1:nrow(patterns), sum(xy$x==x), replace=TRUE, prob=pmat[x,] )
}
xy$y[ xy$sel==1 ] <- sample( 1:nrow(patterns), sum(xy$sel==1), replace=TRUE )

counts <- with(xy, table(x,y) )
nonz <- (counts>0)
nonzcounts <- counts[nonz]

# simple inference
loglik <- function (par) {
    # negative log liklihood
    pmat <- expm( par[1] * singlebase( par[2:(1+nbases*(nbases-1))] ), method="Higham08" )
    return( (-1) * sum( nonzcounts * log(pmat[nonz]) ) )
}

system.time( ans <- optim( par=c(t,baserates), fn=loglik, lower=0, method="L-BFGS-B" ) )
# 22s
# looks good!

# ok, E-M algorithm for nonzero p
em.step <- function (p,pmat) {
    # proportion of each category assigned to "selection" at these parameters
    psel <- p * (1/nrow(patterns)) / ( (1-p)*pmat + p*(1/nrow(patterns)) )
    pp <- ( sum(psel*as.vector(counts))/sum(counts) )
    return(pp)
}

# check:
ppp <- rep(0,30000)
ppp[1] <- runif(1)/100
for (k in 2:length(ppp)) { ppp[k] <- em.step(ppp[k-1],pmat) }
plot(ppp,type='l'); abline(h=p,lty=2)


em.step <- function (par) {
    cat(".")
    pmat <- expm( par[2] * singlebase( par[3:(2+nbases*(nbases-1))] ), method="Higham08" )
    psel <- par[1] * (1/nrow(patterns)) / ( (1-par[1])*pmat + par[1]*(1/nrow(patterns)) )[nonz]
    em.loglik <- function(par) {
        # negative expected log-likelihood
        pmat <- expm( par[2] * singlebase( par[3:(2+nbases*(nbases-1))] ), method="Higham08" )
        ll <- ( (-1) * sum( nonzcounts * ( (1-psel)*log((1-par[1])*pmat[nonz]) + psel*log(par[1]/nrow(patterns)) ) ) )
        if (!is.finite(ll)) { browser() }
        return(ll)
    }
    ans <- optim( par, em.loglik, method="L-BFGS-B", lower=1e-7, upper=c(1-1e-5,rep(Inf,1+nbases*(nbases-1))) )
    return(ans)
}

system.time( em.step(c(p,t,baserates)) )

est.list <- vector("list",40)
est.list[[1]] <- list( par=c( runif(0,.01), runif(0,4), jitter(baserates) ) )
system.time( for (k in seq_along(est.list)[-1]) {
            est.list[[k]] <- em.step( est.list[[k-1]]$par )
        } )

