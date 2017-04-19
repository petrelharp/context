
library(contextual)
library(contextutils)
library(simcontext)

bases <- c("X","O")
mutpats <- list( list(c("X","O"),c("O","X")), list(c("XO","OX")) )
mutrates <- c(1,1)

nsites <- 20
initseq <- rinitseq(nsites,bases)
newseq <- list(simseq(seqlen=nsites,tlen=.05,mutpats=mutpats,mutrates=mutrates,selpats=list(),patlen=2,bases=bases,count.trans=TRUE,initseq=initseq))
for (k in 1:10) { 
    newseq[[k+1]] <- simseq(seqlen=nsites,tlen=.05,mutpats=mutpats,mutrates=mutrates,selpats=list(),patlen=2,bases=bases,count.trans=TRUE,initseq=newseq[[k]]$finalseq)
}
lapply(newseq,show.simseq,latex=TRUE)

nsites <- 20
initseq <- rinitseq(nsites,bases)
newseq <- list(simseq(seqlen=nsites,tlen=.1,mutpats=mutpats,mutrates=mutrates,selpats=list(),patlen=2,bases=bases,count.trans=TRUE,initseq=initseq))
for (k in 1:10) { 
    newseq[[k+1]] <- simseq(seqlen=nsites,tlen=.1,mutpats=mutpats,mutrates=mutrates,selpats=list(),patlen=2,bases=bases,count.trans=TRUE,initseq=newseq[[k]]$finalseq)
}
lapply(newseq,show.simseq,latex=TRUE)

