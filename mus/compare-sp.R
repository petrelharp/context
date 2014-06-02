aligns <- c( "mm9ornAna1", "mm9hg19", "mm9galGal3", "mm9oryCun2", "mm9rn5" )
mle.tables <- lapply( sapply(paste(aligns,"/3.1.counts-results",sep=''),list.files,pattern="*.tsv",full.names=TRUE),
    read.csv, header=TRUE, sep='\t' )
names(mle.tables) <- aligns
mles <- data.frame( t( sapply(mle.tables,function (x) as.matrix(x)["mle",]) ) )
colnames(mles) <- gsub("..","->",gsub("...","|",gsub("tmut..","",colnames(mles)),fixed=TRUE),fixed=TRUE)

mles <- mles[order(mles$branchlen),]
matplot(mles[,1:18],type='l')
legend('topleft',legend=colnames(mles)[1:18],col=1:6,lty=1:5)

