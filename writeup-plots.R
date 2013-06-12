scriptdir <- "./"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))
require(mcmc)


##########
## Tree CpG, with CpG = 0
infile <- "ising/ising-sims/selsims-2013-05-28-17-12-0275615.RData"

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)

lwin <- rwin <- win <- 3
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( "writeup-plots/", basename(basedir), sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

load(datafile)

varnames <- c("lambda","beta","gamma")

all.mrun <- do.call( rbind, lapply( seq_along(mcmcdatafiles), function (k) {
            load(mcmcdatafiles[[k]])
            tmp <- data.frame( mrun$batch )
            names(tmp) <- varnames
            tmp$mcmcnum <- k
            tmp
        } ) )

subseq <- seq(1,nrow(all.mrun),by=max(1,floor(nrow(all.mrun)/1000)))

x <- rbind( all.mrun[subseq,(-1)*match("mcmcnum",colnames(all.mrun))], truth ) 
cols <- c( adjustcolor(rainbow(64),.3)[1+floor(65*(1:(nrow(x)-1))/nrow(x))], adjustcolor("black",.5) )

pdf(file=paste(plotfile,"-traces.pdf",sep=''),width=7,height=5, pointsize=10)
plotthese <- varnames
layout(matrix(1:6,nrow=2))
par(mar=c(4,4,1,1)+.1)
for (k in 1:3) {
    thisone <- match(plotthese[k],names(x))
    thatone <- match(plotthese[1+k%%3],names(x))
    plot( x[,thisone], x[,thatone], xlab=plotthese[k], ylab=plotthese[1+k%%3], col=cols, cex=c(rep(.5,nrow(x)-1),2), pch=20 )
    text(truth[thisone],truth[thatone],"truth", pos=4)
    hist( all.mrun[,thisone], main='', xlab=paste("posterior distribution of",plotthese[k]) )
    abline(v=truth[thisone], col='red')
}
dev.off()

##########
## Tree CpG, with CpG = 0
infile <- "tree-cpg/cpg-tree-sims/selsims-2013-06-03-13-17-0790276.RData"

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)

lwin <- rwin <- 2; win <- 1
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( "writeup-plots/", basename(basedir), sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

load(datafile)

varnames <- c( "branchlen", 
                paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), 
                names(initfreqs)[-length(initfreqs)] 
            )


all.mrun <- do.call( rbind, lapply( seq_along(mcmcdatafiles), function (k) {
            load(mcmcdatafiles[[k]])
            tmp <- data.frame( mrun$batch )
            names(tmp) <- varnames
            tmp$mcmcnum <- k
            tmp
        } ) )

subseq <- seq(1,nrow(all.mrun),by=max(1,floor(nrow(all.mrun)/1000)))

x <- rbind( all.mrun[subseq,(-1)*match("mcmcnum",colnames(all.mrun))], truth ) 
cols <- c( adjustcolor(rainbow(64),.3)[1+floor(65*(1:(nrow(x)-1))/nrow(x))], adjustcolor("black",.5) )

# initial frequencies
pdf(file=paste(plotfile,"-initfreqs.pdf",sep=''),width=7,height=5, pointsize=10)
layout(matrix(1:6,nrow=2))
par(mar=c(4,3,1,1)+.1)
plot( A ~ T, data=x, col=cols, cex=c(rep(.5,nrow(x)-1),2) )
text(.25,.25,"truth", pos=4)
hist( all.mrun$T, main='', xlab="posterior distribution of T" )
abline(v=.25, col='red')
plot( C ~ A, data=x, col=cols, cex=c(rep(.5,nrow(x)-1),2) )
text(.25,.25,"truth", pos=4)
hist( all.mrun$A, main='', xlab="posterior distribution of A" )
abline(v=.25, col='red')
plot( T ~ C, data=x, col=cols, cex=c(rep(.5,nrow(x)-1),2) )
text(.25,.25,"truth", pos=4)
hist( all.mrun$C, main='', xlab="posterior distribution of C" )
abline(v=.25, col='red')
dev.off()

pdf(file=paste(plotfile,"-mutrates.pdf",sep=''),width=7,height=5, pointsize=10)
plotthese <- c( "mut: C->T", "mut: T->C", "mut: CG->TG | CG->CA")
layout(matrix(1:6,nrow=2))
par(mar=c(4,4,1,1)+.1)
for (k in 1:3) {
    thisone <- match(plotthese[k],names(x))
    thatone <- match(plotthese[1+k%%3],names(x))
    plot( x[,thisone], x[,thatone], xlab=plotthese[k], ylab=plotthese[1+k%%3], col=cols, cex=c(rep(.5,nrow(x)-1),2) )
    text(truth[thisone],truth[thatone],"truth", pos=4)
    hist( all.mrun[,thisone], main='', xlab=paste("posterior distribution of",plotthese[k]) )
    abline(v=truth[thisone], col='red')
}
dev.off()
