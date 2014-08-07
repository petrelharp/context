# Test inference from simulation
source("context.R")
source("context-inference-fns.R")
source("sim-context-fns.R")

# do diagnostic plots?
diagnostics <- interactive()

# maximum size of pattern
patlen <- 3
leftwin <- 0
rightwin <- 0

# list of patterns that have mutation rates
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  lapply( combn( bases, 2, simplify=FALSE ), rev),  # single-base rates
        list( c("CG","TG"), c("CG","CA") ), # CpG rates
    NULL )

# mutation rates
mutrates <- runif( length(mutpats) )*1e-8

# list of patterns with selection coefficients
# can be regexps but only using "." and "[...]" (needs regexplen() to work)
selpats <- c(
        "[GC]",
        "[AT]",
        # paste( paste(rep(".",leftwin),collapse=''), paste( codons$codon[codons$aa %in% synons], paste(rep(".",rightwin),collapse=''), sep='' ), sep='' ), # all codons
        # codons$codon[! codons$aa %in% synons],  # stop codons
    NULL )
# selection coefficients
selcoef <- runif( length(selpats) )*1e-4

# check:
stopifnot( patlen >= max( c( sapply(mutpats,function (x) max(nchar(x[1]),nchar(x[2]))), sapply(gsub(".","",selpats,fixed=TRUE),regexplen) ) ) )

# other params
Ne <- 1e4
seqlen <- 1e5
tlen <- 6e6

# simulate the sequence
simseqs <- simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef )


## transition probabilities?
# size of window on either side of the focal site
leftwin <- 2
rightwin <- 2
shortwin <- 1
longwin <- leftwin+shortwin+rightwin

subtransmatrix <- gettransmatrix(mutpats*tlen, mutrates, selpats, selcoef, Ne, shortwin, leftwin, rightwin)
counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs=simseqs, leftwin=leftwin )

# averaged
expected <- (seqlen-longwin+1) * (1/nbases)^longwin * subtransmatrix
# accounting for initial sequence
in.expected <- rowSums(counts) * subtransmatrix
# indicator of counts where patterns have changed
changed <- whichchanged(subtransmatrix,leftwin=leftwin,shortwin=shortwin)

# Compare observed and expected counts:
layout(matrix(1:4,2))
plot( as.vector(counts), as.vector(expected), col=1+changed ); abline(0,1)
plot( as.vector(counts), as.vector(in.expected), col=1+changed ); abline(0,1)
plot( as.vector(counts)[changed], as.vector(expected)[changed], col=2 ); abline(0,1)
plot( as.vector(counts)[changed], as.vector(in.expected)[changed], col=2 ); abline(0,1)
# should be mixture of poissons -- compare means in bins
usethese <- ( changed & as.vector(in.expected)>.1 )
exp.bins <- cut( as.vector(in.expected)[usethese], 40 )
exp.bin.vals <- tapply( as.vector(in.expected)[usethese], exp.bins, mean )
points( tapply( as.vector(counts)[usethese], exp.bins, mean ), exp.bin.vals, pch=20 )

# how many events per window?
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=4*longwin))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=2*longwin))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=longwin))
hist(simseqs$ntrans$loc,breaks=seq(0,seqlen+10,by=1))



####
# Inference?

leftwin <- 2
rightwin <- 2
shortwin <- 1
longwin <- leftwin+shortwin+rightwin
genmatrix <- makegenmatrix( mutpats, selpats, patlen=longwin)
genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )

subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs, leftwin )

nmuts <- length(mutpats)
nsel <- length(selpats)
const <- sum(counts*log(1/1.09))  # rescale to get value closer to zero
rate.scale <- max(rowSums(genmatrix))
sel.scale <- 1e-5
Ne.scale <- 1e4
likfun <- getlikfun(nmuts=nmuts,nsel=nsel,genmatrix=genmatrix,projmatrix=projmatrix,const=const)

initpar <- c( 
        runif( length(mutpats) )*1e-8 * tlen, 
        runif( length(selpats) )*1e-4 / coef.scale,
        runif(1)*2*Ne / Ne.scale ,
        ) # falsehood
initpar <- c( runif( length(mutpats) ) * mean(mutrates) * tlen, runif( length(selpats) )*1e-4, runif(1)*2*Ne ) # random init
truth <- c( mutrates * tlen, selcoef, Ne )  # truth
lbs <- c( rep( 1e-6, nmuts+nsel ), 1e-2 )
ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs )
lbs <- c( rep( 1e-6, nmuts+nsel ), 1e-2 )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=c( rep(1,nmuts), rep(sel.scale,nsel), Ne.scale ) ) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=c( rep(1,nmuts), rep(sel.scale,nsel), Ne.scale ) ) )

estimates <- cbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth )
apply( estimates, 2, likfun )

#####
# Now, evolve from a common root.

seqlen <- 1e5
tlen <- 6e6
Ne <- 1e4
initseq <- rinitseq(seqlen,bases)
simseqs <- lapply(1:2, function (k) simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef ) )

## transition probabilities?
# size of window on either side of the focal site
leftwin <- 2
rightwin <- 2
shortwin <- 1
longwin <- leftwin+shortwin+rightwin

subtransmatrix <- gettransmatrix(mutpats, mutrates*tlen, selpats, selcoef, Ne, shortwin, leftwin, rightwin)
for (k in seq_along(simseqs)) {
    simseqs[[k]]$counts <- counttrans( rownames(subtransmatrix), colnames(subtransmatrix), simseqs=simseqs[[k]], leftwin=leftwin )
}


#######
# testing: number of subsequences

xseqlen <- 1000
xpats <- getpatterns(2)
x <- replicate( 1000, { 
            aseq <- paste( sample(bases,xseqlen,replace=TRUE), collapse="" )
            xmatches <- sapply( lapply( xpats, function (p) gregexpr(paste("(?=",p,")",sep=''),aseq,perl=TRUE) ), "[[", 1 )
            sapply( xmatches, length )
        } )
rownames(x) <- xpats
poisx <- rpois(10000,lambda=xseqlen/length(xpats))
hx <- hist(c(x,poisx), freq=FALSE, plot=FALSE)
histx <- invisible( lapply( seq_along(xpats), function (k) hist( x[k,], col=adjustcolor(rainbow(64)[k],.2), plot=FALSE, breaks=hx$breaks, freq=FALSE ) ) )
matplot( sapply(histx, "[[", "density" ), type='l', lty=c(1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,1), col=1:16 )
legend("topright",legend=xpats,lty=c(1,2,2,2,2,1,2,2,2,2,1,2,2,2,2,1),col=1:16)
lines( hist( poisx, breaks=hx$breaks, plot=FALSE )$density, lwd=2 )
