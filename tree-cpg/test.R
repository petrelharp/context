f <- function(params) {
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:nfreqs)]
    pred <- predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[1]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=rev(branchlens) )
    (counts[[1]]-pred)
}

do.call(cbind,lapply(mutpats[-13],"[[",1))

dpAT <- 1e-3*c( rep(0,14), c(1,-1), c(0,0) )
dmAT <- rep(0,18); dmAT[1+c(1,7)] <- 1e-4*c(1,-1)  # A->T, T->A
dzAT <- rep(0,18); dzAT[1+c(1,7)] <- 1e-3*c(1,-1)*params[1+c(1,7)]/.25  # A->T, T->A

dpAT <- 1e-3*c( rep(0,14), c(1,-1), c(0,0) )
dmAT <- rep(0,18); dmAT[1+c(1,7)] <- 1.5e-4*c(1,-1)  # A->T, T->A
dzAT <- 1.5e-5*c(0,1,1,1,-1,-1,0,-1,0,0,0,0,0,0,0,0,0)

c( dp=lud(params+dpAT), dm=lud(params+dpAT+dmAT), dz=lud(params+dpAT+dmAT+dzAT) ) - lud(params)

orig <- f(params)
modp <- f(params+dpAT)
modm <- f(params+dpAT+dmAT)
modz <- f(params+dpAT+dmAT+dzAT)
labels <- c("AATAA-T","AAAAA-T","AACAA-T","GGGGG-A","GGTGG-A","GGAGG-A","CCACC-C","GAATC-G","TTTAC-C","CTTTC-G")
alllabs <- outer(rownames(f(params)),colnames(f(params)),paste,sep="-")

plot(modp-orig,ylim=range(as.vector(modp-orig),as.vector(modm-orig),as.vector(modz-orig)))
points(modm-orig,col=adjustcolor('green',.25))
points(modz-orig,col=adjustcolor('red',.25))
text( match(labels,alllabs), (modp-orig)[match(labels,alllabs)], labels=labels, pos=c(2,2,2,4,4,4,3,3,1,1) )
text( match(labels,alllabs), (modm-orig)[match(labels,alllabs)], labels=labels, pos=c(2,2,2,4,4,4,3,3,1,1), col='green' )
text( match(labels,alllabs), (modz-orig)[match(labels,alllabs)], labels=labels, pos=c(2,2,2,4,4,4,3,3,1,1), col='red' )

identify(seq_along(f(params)),f(params+dpAT)-f(params),alllabs)

########
tmp <- lud(params)
c( mut=lud(params+dmAT)-tmp,
    pi=lud(params+dpAT)-tmp,
    both=lud(params+dpAT+dmAT)-tmp  
    )

tmp <- lud(truth)
c( mut=lud(truth+dmAT)-tmp,
    pi=lud(truth+dpAT)-tmp,
    both=lud(truth+dpAT+dmAT)-tmp  
    )
