library(expm)
nstates <- 16
ndiag <- 10

M <- matrix( 0, nstates, nstates )
for (k in 1:ndiag) {
    M[ row(M)==col(M)+k ] <- runif(nstates-k)
    M[ row(M)==col(M)-k ] <- runif(nstates-k)
}
diag(M) <- -rowSums(M)
tt <- 1

# Pt <- lapply(seq_along(tt), function (k) expm(diff(c(0,tt))[k]*M) )
Pt <- expm(tt*M)

# data
xx <- data.frame(x0=sample(1:nstates,1000000,replace=TRUE))
for (k in seq_along(tt)) {
#     xx <- cbind(xx, sapply(xx[,k], function (x) sample(1:nstates,size=1,prob=Pt[[k]][x,])) )
     xx <- cbind(xx, sapply(xx[,k], function (x) sample(1:nstates,size=1,prob=Pt[x,])) )
}
names(xx) <- paste("x",0:length(tt),sep='')

obs <- table(xx$x0,xx$x1)
nonzobs <- obs!=0
rm(xx)

########
# inference with expm
nonzM <- (abs(row(M)-col(M)) <= ndiag ) & (row(M) != col(M))

fn <- function (z) {
    zz <- matrix(0,nstates,nstates)
    zz[nonzM] <- z
    diag(zz) <- -rowSums(zz)
    fnz <- (-1) * sum( obs[nonzobs] * log( expm(zz)[nonzobs] ) ) 
    if (!is.finite(fnz)) { recover() } else { return(fnz) }
}

initpar <- tt*M[nonzM]
system.time( ans <- optim( par=initpar, fn=fn, lower=0, method="L-BFGS-B" ) )
# 28 sec for nstates=16
# 106 sec for nstates=32
# 337 sec for nstates=48
# 669 sec for nstates=64
# predicts (.4361 * nstates - 2.5)^2 seconds
# predicts 47 min for nstates=128
# predicts 198 min for nstates=256
Rprof("expmprof.txt")
ans <- optim( par=initpar, fn=fn, lower=0, method="L-BFGS-B" )
Rprof(NULL)
summaryRprof("expmprof.txt")

estMt <- matrix(0,nstates,nstates)
estMt[nonzM] <- ans$par
diag(estMt) <- (-1)*rowSums(estMt)
estPt <- expm(estMt)

layout(matrix(1:4,2,2))
plot( as.vector(Pt), as.vector(estPt), main ="Pt" )
abline(0,1)
matplot( Pt, type='l', col=adjustcolor("black",.7) )
matlines( estPt, col=adjustcolor("red",.7) )
plot( as.vector(tt*M)[nonzM], as.vector(estMt)[nonzM], main="t*M" )
abline(0,1)
matplot( tt*M-diag(diag(tt*M)), type='l', col=adjustcolor("black",.7), log='y' )
matlines( estMt-diag(diag(estMt)), col=adjustcolor("red",.7) )

########
# inference without expm
#  deal with row sums by putting penalty on row sum minus 1

w <- obs[nonzobs]/mean(obs[nonzobs])
obsrow <- rowSums(obs)/mean(obs[nonzobs])
nonzM <- (row(M) != col(M))
fn <- function (z) {
    dim(z) <- c(nstates,nstates)
    rsz <- rowSums(z)
    fnz <- -sum( w * log( z[nonzobs] ) ) + sum(obsrow*log(rsz))  + sum( (rsz-1)^2 )
    if (!is.finite(fnz)) { recover() } else { return(fnz) }
}

initpar <- Pt
system.time( ans <- optim( par=initpar, fn=fn, lower=1e-6, method="L-BFGS-B" ) )
dim(ans$par) <- c(nstates,nstates)
estPt <- sweep(ans$par,1,rowSums(ans$par),"/")

# 64 states = 205 seconds or 114 seconds

layout(1:2)
plot( as.vector(Pt), as.vector(estPt) )
abline(0,1)
matplot( Pt, type='l', col=adjustcolor("black",.5) )
matlines( estPt, col=adjustcolor("red",.5) )
matlines( sweep(obs,1,rowSums(obs),"/"), col=adjustcolor("green",.5) )  # this is the correct MLE

######
# triad inference

nstates <- 32
ndiag <- 10

M <- matrix( 0, nstates, nstates )
for (k in 1:ndiag) {
    M[ row(M)==col(M)+k ] <- runif(nstates-k)
    M[ row(M)==col(M)-k ] <- runif(nstates-k)
}
diag(M) <- -rowSums(M)
tt <- c(1,1.5,2)

Pt <- lapply(seq_along(tt), function (k) expm(tt[k]*M) )

# data
initprob <- rep(1/nstates,nstates)
xx <- data.frame(x0=sample(1:nstates,1000000,replace=TRUE,prob=initprob))
for (k in seq_along(tt)) {
     xx <- cbind(xx, sapply(xx[,k], function (x) sample(1:nstates,size=1,prob=Pt[[k]][x,])) )
}
names(xx) <- paste("x",0:length(tt),sep='')

obs <- table(xx[2:length(xx)])
nonzobs <- obs!=0
rm(xx)

nonzM <- (abs(row(M)-col(M)) <= ndiag ) & (row(M) != col(M))
w <- obs[nonzobs]/mean(obs[nonzobs])

fn <- function (z) {
    # z is ( tt, initprob, M[nonzM] )
    tt <- z[1:length(tt)]
    ip <- z[(length(tt)+1):(length(tt)+nstates)]
    zz <- matrix(0,nstates,nstates)
    zz[nonzM] <- z[(length(tt)+nstates+1):length(z)]
    dim(zz) <- c(nstates,nstates)
    diag(zz) <- -rowSums(zz)
    PP <- lapply(tt, function (t) expm(t*zz))
    P <- colSums( ip * PP[[1]][,rep(1:nstates,each=nstates^2)] * PP[[2]][,rep(1:nstates,nstates,each=nstates)] * PP[[3]][,rep(1:nstates,nstates^2)] )
    fnz <- (-1) * sum( w * log( P[nonzobs] ) ) + (sum(ip)-1)^2
    if (!is.finite(fnz)) { recover() } else { return(fnz) }
}

initpar <- c(tt,initprob,M[nonzM])
system.time( ans <- optim( par=initpar, fn=fn, lower=0, method="L-BFGS-B" ) )
# 51 sec for nstates=16

esttt <- ans$par[seq_along(tt)]
estip <- ans$par[(length(tt)+1):(length(tt)+nstates)]
estip <- estip/sum(estip)
estMt <- matrix(0,nstates,nstates)
estMt[nonzM] <- ans$par[(length(tt)+nstates+1):length(ans$par)]
diag(estMt) <- (-1)*rowSums(estMt)

layout(matrix(1:4,2,2))
plot( c(tt,nstates*initprob), c(esttt,nstates*estip), col=c(rep("red",length(tt)),rep('black',length(estip))) )
abline(0,1)
plot( as.vector(M), as.vector(estMt) )
abline(0,1)
plot( as.vector(M[nonzM]), as.vector(estMt[nonzM]) )
abline(0,1)
matplot( M, type='l', col=adjustcolor("black",.5) )
matlines( estMt, col=adjustcolor("red",.5) )

