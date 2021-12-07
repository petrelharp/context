chrs <- read.table("hu_ch_gor_differences_chr1.counts.gz",header=TRUE,stringsAsFactors=TRUE)
names(chrs)[4] <- "chr1"
for (k in 2:22) { 
    chrs <- merge(chrs,read.table(paste("hu_ch_gor_differences_chr",k,".counts.gz",sep=''),header=TRUE,stringsAsFactors=TRUE),by=c("Hu_Quintet","Ch_Quintet","Gor_Quintet"),all=FALSE) 
    names(chrs)[match("count",names(chrs))] <- paste("chr",k,sep='') 
}
chrs <- subset(chrs, (!grepl("[^ACGT]",Hu_Quintet)) & (!grepl("[^ACGT]",Ch_Quintet)) & (!grepl("[^ACGT]",Gor_Quintet)) )
for (k in 1:3) { chrs[,k] <- droplevels(chrs[,k]) }
save(chrs,file="all-chr-counts.RData")

allpats <- sort(unique(do.call(c,lapply(1:3,function(k)levels(chrs[,k])))))
pats <- chrs[,1:3]
for (k in 1:3) { pats[,k] <- factor( levels(chrs[,k])[chrs[,k]], levels=allpats ) } 
counts <- chrs[,4:25]
y <- sweep(counts,2,colSums(x),"/")
unchanged <- ( (pats[,1] == pats[,2]) & (pats[,2] == pats[,3]) )
gcpats <- sapply( strsplit( allpats, "" ), function (x) sum( x %in% c("C","G") ) )
gccontent <- with( pats, gcpats[Hu_Quintet] + gcpats[Ch_Quintet] + gcpats[Gor_Quintet] )/3
gccols <- rainbow(32)[cut(gccontent,24)]
count.changes <- function (...) { 
    # count number of changed sites
    x <- list(...)
    xx <- lapply(x, function (y) {
            yy <- unlist( sapply( as.character(y), strsplit, '') )
            dim(yy) <- c(length(yy)/length(y), length(y))
            yy
        } )
    zz <- ( xx[[1]]!=xx[[2]] )
    for (k in seq_along(xx[-(1:2)])) {
        zz <- ( zz | ( xx[[2+k]] != xx[[1]] ) )
    }
    colSums(zz)
}
nchanges <- with( pats, data.frame(
            all=count.changes(Hu_Quintet,Ch_Quintet,Gor_Quintet) ,
            hu_ch=count.changes(Hu_Quintet,Ch_Quintet) ,
            hu_go=count.changes(Hu_Quintet,Gor_Quintet) ,
            ch_go=count.changes(Ch_Quintet,Gor_Quintet) 
        ) )

matplot(y,type='l',log='y',col=adjustcolor(rainbow(22),.2))

# correlations between chromosomes and GC content
lmat <- matrix( 1:22^2, nrow=22 )
lmat[upper.tri(lmat,diag=FALSE)] <- 1:choose(22,2)
lmat[lower.tri(lmat,diag=TRUE)] <- 0

png(file="all-vs-all-unchanged.png",width=30*144,height=30*144,res=144)
layout( lmat )
par(mar=c(0,0,0,0))
for (k in 1:21) { for (j in (k+1):22) { 
    plot( y[unchanged,k], y[unchanged,j], xaxt='n', yaxt='n', xlab='', ylab='', pch=20, col=gccols[unchanged] ) 
    abline(0,1,lwd=2)
} }
dev.off()

png(file="all-vs-all-changed.png",width=30*144,height=30*144,res=144)
layout( lmat )
par(mar=c(0,0,0,0))
for (k in 1:21) { for (j in (k+1):22) { plot( chrs[!unchanged,3+k], chrs[!unchanged,3+j], xaxt='n', yaxt='n', xlab='', ylab='', pch=".", col=gccols[!unchanged] ) } }
dev.off()

# number of changes
tapply( rowSums(counts), nchanges$all, sum )
tapply( rowSums(counts), nchanges[2:4], sum, na.rm=TRUE )
tapply( rowSums(counts)/sum(counts), nchanges[2:4], sum, na.rm=TRUE )

