# Test inference from simulation
source("codons.R")
source("codon-inference-fns.R")
source("sim-context-fns.R")

# maximum size of pattern (for simulation)
patlen <- 2
lwin <- 0
rwin <- 0
mutpats <- c( 
        combn( bases, 2, simplify=FALSE ), 
        lapply( combn( bases, 2, simplify=FALSE ), rev),
        list(
                gc=list( c("GC","TC") ) ,
                ga=list( c("GA","CA") ) ,
                ta=list( c("AT","GT") ) 
            )
    ) 
mutrates <- runif( length(mutpats) )*3e-8
selpats <- c( "[GC]", "[AT]", "AAA", "CCC" )
selcoef <- sort( runif( length(selpats) ) )*1e-4
# other params
Ne <- 1e4
seqlen <- 1e4
tlen <- 6e7  # 6e7 gives lots of transitions; 1e7 not so many
# simulate the sequence
simseqs <- simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef )

####
# Inference.
lwin <- 0
rwin <- 2
win <- 1
winlen <- lwin+win+rwin

genmatrix <- makegenmatrix( mutpats, selpats, patlen=winlen)
genmatrix@x <- update(genmatrix,mutrates*tlen,selcoef,Ne)
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )

subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs, lwin )

nmuts <- length(mutpats)
nsel <- length(selpats)
const <- sum(counts*log(1/1.09))  # rescale to get value closer to zero
rate.scale <- max(rowSums(genmatrix))
sel.scale <- 1e-5
Ne.scale <- 1e4
# params are: mutrates, selcoef, Ne, *scaled*
likfun <- getlikfun(nmuts=nmuts,nsel=nsel,genmatrix=genmatrix,projmatrix=projmatrix,const=const)
initpar <- c( runif( length(mutpats) ) * mean(mutrates) * tlen, runif( length(selpats) )*1e-4, runif(1)*2*Ne ) # random init
truth <- c( mutrates * tlen, selcoef, Ne )  # truth
lbs <- c( rep( 1e-8, nmuts+nsel ), 1e2 )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=c( rep(1,nmuts), rep(sel.scale,nsel), Ne.scale ) ) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=c( rep(1,nmuts), rep(sel.scale,nsel), Ne.scale ) ) )

estimates <- rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth )
estimates
apply( estimates, 1, likfun )

# look at param estimates
layout(1:2)
matplot(t((estimates-estimates[,"truth"])/rowMeans(estimates)), type='l',col=1:4,lty=1:4)
abline(v=.5+c(nmuts,nmuts+nsel))
legend("topleft",col=1:4,lty=1,legend=rownames(estimates))
matplot(t(estimates/truth),type='l',col=1:4,lty=1:4)
abline(v=.5+c(nmuts,nmuts+nsel))
legend("topleft",col=1:4,lty=1,legend=rownames(estimates))

# look at observed and expected.
true.expected <- rowSums(counts) * subtransmatrix
est.genmatrix <- genmatrix
est.genmatrix@x <- update( genmatrix, mutrates=estimates["ans",1:nmuts], selcoef=estimates["ans",nmuts+(1:nsel)], Ne=estimates["ans",nmuts+nsel+1] )
est.expected <- rowSums(counts) * computetransmatrix( est.genmatrix, projmatrix=projmatrix )
cheating.genmatrix <- genmatrix
cheating.genmatrix@x <- update( genmatrix, mutrates=estimates["cheating",1:nmuts], selcoef=estimates["cheating",nmuts+(1:nsel)], Ne=estimates["cheating",nmuts+nsel+1] )
cheating.expected <- rowSums(counts) * computetransmatrix( cheating.genmatrix, projmatrix=projmatrix )

layout(matrix(1:4,nrow=2))
for (k in 1:4) {
    lord <- order( true.expected[,k] )
    plot( counts[lord,k], xaxt='n', xlab='', main=colnames(counts)[k] )
    axis(1,at=1:nrow(counts),labels=rownames(counts)[lord],las=3)
    lines(true.expected[lord,k],col='red')
    lines(est.expected[lord,k],col='green')
    lines(cheating.expected[lord,k],col='blue',lty=2,lwd=2)
    legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
}

# indicator of counts where patterns have changed
changed <- whichchanged(subtransmatrix,lwin=lwin,win=win)
