scriptdir <- "./"
source(paste(scriptdir,"context-inference-fns.R",sep=''),chdir=TRUE)
source(paste(scriptdir,"sim-context-fns.R",sep=''),chdir=TRUE)
require(mcmc)


##########
## Ising
infile <- "ising/ising-sims/selsims-2013-05-28-17-12-0275615.RData"

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
# translate from old version
lwin <- leftwin; rwin <- rightwin; win <- shortwin; winlen <- longwin

leftwin <- rightwin <- shortwin <- 3
basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( "writeup-plots/", basename(basedir), sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

load(datafile)
# translate from old version
lwin <- leftwin; rwin <- rightwin; win <- shortwin; winlen <- longwin

varnames <- c("lambda","beta","gamma")

all.mrun <- do.call( rbind, lapply( seq_along(mcmcdatafiles), function (k) {
            load(mcmcdatafiles[[k]])
            tmp <- data.frame( mrun$batch )
            names(tmp) <- varnames
            tmp$mcmcnum <- k
            tmp
        } ) )

win.1.1 <- which( grepl("1-1-1",mcmcdatafiles) & table( all.mrun$mcmcnum ) >= 10000 )
win.3.0 <- which( grepl("0-3-0",mcmcdatafiles) & table( all.mrun$mcmcnum ) >= 10000 )
win.3.3 <- which( grepl("3-3-3",mcmcdatafiles) & table( all.mrun$mcmcnum ) >= 10000 )

subseq <- seq( min(which(all.mrun$mcmcnum %in% win.3.3)), max(which(all.mrun$mcmcnum %in% win.3.3)), by=max(1,floor(sum(all.mrun$mcmcnum %in% win.3.3)/1000)) )

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
    hist( all.mrun[all.mrun$mcmcnum %in% win.3.3,thisone], main='', xlab=paste("posterior distribution of",plotthese[k]) )
    abline(v=truth[thisone], col='red')
}
dev.off()

pdf(file=paste(plotfile,"-estimate-hists.pdf",sep=''),width=3,height=3,pointsize=10)
par(mar=c(4,3,.3,0)+.1,mgp=c(2.1,1,0))
layout(t(1:3))
for (mcnums in list( win.3.3, win.3.0, win.1.1 ) ) {
tmp <- sweep(abs(all.mrun)[all.mrun$mcmcnum %in% mcnums,1:3],2,truth,"/")
names(tmp) <- varnames
hist(tmp[,1],breaks=50,col=adjustcolor('grey',.5), border=adjustcolor("black",.5), main='', xlim=c(.9,1.1), ylab='posterior density', freq=FALSE, xlab="" )
hist(tmp[,2],breaks=50,col=adjustcolor("blue",.5),border=adjustcolor("black",.25), add=TRUE, freq=FALSE)
hist(tmp[,3],breaks=50,col=adjustcolor("red",.5),border=adjustcolor("black",.5), add=TRUE, freq=FALSE)
abline(v=1,lwd=2)
legend( "topleft", fill=c(adjustcolor('grey',.5),adjustcolor("blue",.5),adjustcolor("red",.5)), legend=as.expression(c(substitute(lambda),substitute(beta),substitute(gamma))) )
}
dev.off()

# for talk
pdf(file=paste(plotfile,"-estimate-boxplots.pdf",sep=''),width=3,height=3,pointsize=10)
par(mar=c(5,3,.3,0)+.1,mgp=c(2.1,1,0))
tmp <- with( subset(all.mrun,mcmcnum%in%c(win.3.3, win.3.0, win.1.1)), 
    data.frame( val=abs(c(10*lambda,beta,gamma)), mcmcnum=rep(mcmcnum,3), var=rep(c("lambda","beta","gamma"),each=length(lambda)) ) )
tmp$window <- ""; 
tmp$window[tmp$mcmcnum %in% win.3.3] <- "3-3-3"
tmp$window[tmp$mcmcnum %in% win.1.1] <- "1-1-1"
tmp$window[tmp$mcmcnum %in% win.3.0] <- "0-3-0"
tmp$window <- factor(tmp$window,levels=c("0-3-0","1-1-1","3-3-3"))
boxplot( val ~ var + window, at=c(1:3,5:7,9:11), data=tmp, las=3, col=adjustcolor(c("black","blue","red"),.5), xaxt='n', ylim=c(.85,1.1) )
axis(1, at=c(2,6,9), labels=levels(tmp$window), tick=FALSE)
abline(h=1,col='red',lwd=2)
abline(v=c(4,8),lty=2,col='grey')
legend( "topleft", fill=c(adjustcolor('grey',.5),adjustcolor("blue",.5),adjustcolor("red",.5)), legend=as.expression(c(substitute(lambda),substitute(beta),substitute(gamma))), bg='white' )
dev.off()


##########
## Tree CpG, with CpG = 0
infile <- "tree-cpg/cpg-tree-sims/selsims-2013-06-03-13-17-0790276.RData"

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
# translate from old version
lwin <- leftwin; rwin <- rightwin; win <- shortwin; winlen <- longwin

leftwin <- rightwin <- 2; shortwin <- 1
basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( "writeup-plots/", basename(basedir), sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

load(datafile)
# translate from old version
lwin <- leftwin; rwin <- rightwin; win <- shortwin; winlen <- longwin

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

##########
## Tree CpG with CpG larger
infile <- "tree-cpg/cpg-tree-sims/selsims-2013-06-03-13-17-0187525.RData"

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
# translate from old version
lwin <- leftwin; rwin <- rightwin; win <- shortwin; winlen <- longwin

leftwin <- rightwin <- 2; shortwin <- 1
basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( "writeup-plots/", basename(basedir), sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

load(datafile)
# translate from old version
lwin <- leftwin; rwin <- rightwin; win <- shortwin; winlen <- longwin

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


pdf(file=paste(plotfile,"-estimate-boxplots.pdf",sep=''),width=5,height=3,pointsize=10)
layout(t(1:2))
par(mar=c(6,3,.3,0)+.1,mgp=c(2.1,1,0))
tmp <- subset(all.mrun,mcmcnum==5)
names(tmp) <- gsub("mut: ","",names(tmp))
names(tmp)[1+13] <- "CG->TG/CA"
names(tmp)[1] <- "branch len"
boxplot( tmp[,2:13], las=3, ylab="scaled mutation rate", ylim=c(0,.2) )
abline(h=.15,col='red',lwd=2)
boxplot( tmp[,c(1+13,1,1+14:16)] ,las=3, ylim=c(.2,.5))
segments(x0=(1:5)-.5,x1=(1:5)+.5,y0=truth[c(1+13,1,1+14:16)],col='red',lwd=2)
dev.off()

