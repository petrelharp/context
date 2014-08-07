# Evolve from a common root
source("context.R")
source("context-inference-fns.R")
source("sim-context-fns.R")

# maximum size of pattern (for simulation)
patlen <- 2
mutpats <- list( 
        fwd=combn( bases, 2, simplify=FALSE ), 
        rev=lapply( combn( bases, 2, simplify=FALSE ), rev),
        gc=list( c("GC","TC") ) ,
        ga=list( c("GA","CA") ) ,
        ta=list( c("AT","GT") ) 
    ) 
mutrates <- runif( length(mutpats) )*1e-8
selpats <- c( "[GC]", "[AT]" )
selcoef <- runif( length(selpats) )*1e-3
# other params
Ne <- c(1e4,1e4)
seqlen <- 1e4
tlen <- 3e7  # total length 6e7 gives lots of transitions (but still signal); 1e7 not so many
branchlens <- c(1,1)

# tlen <- .1e7: 2.9 sec ; 3e7: 188 sec

initfreqs <- c(.3,.3,.2,.2)
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- lapply(1:2, function (k) simseq( seqlen, tlen*branchlens[k], patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, Ne=Ne[k], initseq=initseq ) ) 
    )

thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- format(Sys.time(), "%Y-%m-%d-%H-%M")
save( thisone, now, patlen, mutpats, selpats, selcoef, Ne, seqlen, tlen, branchlens, initfreqs, simseqs, file=paste(now,thisone,"selsims.RData",sep='') )


####
# Inference.
leftwin <- 0
rightwin <- 2
shortwin <- 1
longwin <- leftwin+shortwin+rightwin

genmatrix <- makegenmatrix( mutpats, selpats, patlen=longwin)
genmatrix@x <- update(genmatrix,mutrates*tlen,selcoef,Ne)
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$finalseq, simseqs[[2]]$finalseq, leftwin=leftwin ),
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[2]]$finalseq, simseqs[[1]]$finalseq, leftwin=leftwin )
        )

nmuts <- length(mutpats)
nsel <- length(selpats)
sel.scale <- 1e-5
Ne.scale <- 1e4
# params are: mutrates, selcoef, Ne, *scaled*

# move from base frequencies (what we estimate) to pattern frequencies
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )
onevec <- rep(1,longwin)
likfun <- function (params) {
        # params are: Ne, branchlens[-1], mutrates, selcoef, basefreqs
        #  -> assume branches are equal length for now
        Ne <- params[1:2]
        branchlens <- c(1,params[3])
        mutrates <- params[4+(1:nmuts)]
        selcoef <- params[4+nmuts+(1:nsel)]
        initfreqs <- params[4+nmuts+nsel+(1:4)][ patcomp ] 
        dim(initfreqs) <- dim(patcomp)
        initfreqs <- as.vector( initfreqs %*% onevec )
        # this is collapsed transition matrix
        updownbranch <- getupdowntrans( genmatrix, projmatrix, list(mutrates,mutrates), list(selcoef,selcoef), c(Ne,Ne), initfreqs, tlens=branchlens )
        # return negative log-likelihood 
        (-1) * ( sum( counts[[1]] * log(updownbranch) ) + sum( counts[[2]] * log(updownbranch) ) )
}

initpar <- c( runif(1)*2*Ne, runif(2)*2*branchlens, runif( length(mutpats) ) * mean(mutrates) * tlen, runif( length(selpats) )*1e-4, rep(1/4,4) ) # random init
truth <- c( Ne,  branchlens, mutrates * tlen, selcoef, initfreqs )  # truth
lbs <- c( c(1,2), rep(1e-3,2), rep( 1e-6, nmuts+nsel), rep(0,4) )
scalevec <- c( rep(Ne.scale,2), rep(1,2), rep(1,nmuts), rep(sel.scale,nsel), rep(1,4) )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=scalevec) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=scalevec) )

estimates <- rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth )
estimates
apply( estimates, 1, likfun )

# look at observed and expected.
updownbranch <- getupdowntrans( genmatrix, projmatrix, list(mutrates*tlen,mutrates*tlen), list(selcoef,selcoef), c(Ne,Ne), initfreqs=rep(1/nrow(genmatrix),nrow(genmatrix)), tlens=c(1,1) )
true.expected <- lapply( counts, function (cnt) { rowSums(cnt) * updownbranch } )
est.updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=list(estimates["ans",1:nmuts],estimates["ans",1:nmuts]), 
    selcoef=list(estimates["ans",nmuts+(1:nsel)],estimates["ans",nmuts+(1:nsel)]), 
    Ne=c(estimates["ans",nmuts+nsel+1],estimates["ans",nmuts+nsel+1]),
    initfreqs=rep(1/nrow(genmatrix),nrow(genmatrix)), tlens=c(1,1) )
est.expected <- lapply( counts, function (cnt) { rowSums(cnt) * est.updownbranch } )
cheating.updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=list(estimates["cheating",1:nmuts],estimates["cheating",1:nmuts]), 
    selcoef=list(estimates["cheating",nmuts+(1:nsel)],estimates["cheating",nmuts+(1:nsel)]), 
    Ne=c(estimates["cheating",nmuts+nsel+1],estimates["cheating",nmuts+nsel+1]),
    initfreqs=rep(1/nrow(genmatrix),nrow(genmatrix)), tlens=c(1,1) )
cheating.expected <- lapply( counts, function (cnt) { rowSums(cnt) * cheating.updownbranch } )
init.counts <- rowSums(counttrans(rownames(genmatrix),colnames(projmatrix),simseqs= simseqs[[1]]))

layout(matrix(1:4,nrow=2))
for (k in 1:4) {
    lord <- order( true.expected[[1]][,k] + true.expected[[2]][,k] )
    plot( counts[[1]][lord,k], type='n', xaxt='n', xlab='', main=colnames(counts[[1]])[k], ylim=range(c(unlist(true.expected),unlist(lapply(counts,as.matrix)),unlist(est.expected))) )
    axis(1,at=1:nrow(counts[[1]]),labels=rownames(counts[[1]])[lord],las=3)
    for (j in 1:2) {
        points( counts[[j]][lord,k], pch=j )
        lines(true.expected[[j]][lord,k],col='red', lty=j)
        lines(est.expected[[j]][lord,k],col='green', lty=j, lwd=2)
        lines(cheating.expected[[j]][lord,k],col='blue',lty=j)
    }
    legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
}

