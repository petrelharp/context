library(contextual)

mcmcs <- list.files("RegulatoryFeature-regions-from-axt/", 
            "biochem-v3-fit-mcmc.*RData", recursive=TRUE, full.names=TRUE)

source("models/mutrate_names.R") # provides mutrate.names

mcmcs <- mcmcs[sapply(mcmcs, function (mcfile) { load(mcfile); (dim(model@results[["batch"]])[1]>1000) } )]
load(mcmcs[1])
mnames <- names(mutrate.names)[match(gsub(".*[.]","",names(coef(model))[model@results$use.par]), mutrate.names)]

png(file="models/biochem-v3.mcmc_traces/biochem-v3.mcmc_traces_%d.png", 
    width=11*300, height=8*300, res=300, type='cairo')
for (mcfile in mcmcs) {
    load(mcfile)
    {
        layout(matrix(c(1,2,3,3), nrow=2), width=c(4,1))
        par(mar=c(0,3,3,1)+.1)
        matplot(model@results[["batch"]], type='l', 
            xaxt='n', xlab='', main=gsub("R[^/]*/*","",dirname(dirname(mcfile))),
            ylab='parameter value')
        par(mar=c(4,3,0,1)+.1)
        matplot(model@results[["batch"]], type='l', ylim=c(-1,2),
            xlab='MCMC step', ylab='parameter value')
        par(mar=c(4,0,3,1)+.1)
        plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n')
        legend("topright", lty=1:5, col=1:6, legend=mnames)
    }
}
dev.off()

xspace <- 4
boxes <- do.call(cbind, lapply(lapply(mcmcs, function (mcfile) { load(mcfile); boxplot(as.data.frame(model@results[["batch"]][seq(20001,40000),]), plot=FALSE) } ), "[[", "stats"))
atvecs <- unlist(lapply(seq_along(mcmcs), function (k) { (seq_along(mnames)-1)*(8+xspace) + k }))
colvec <- rep(1:8, each=length(mnames))
mcnames <- gsub("noOverlap", "nongene", gsub("overlap", "gene", gsub("-knownGeneTx", "", sapply(mcmcs, function (mcfile) gsub("/", " : ", gsub("R[^/]*[/]*","",dirname(dirname(mcfile))) )))))

typecols <- RColorBrewer::brewer.pal("Paired", n=8)[c(1,3,5,7,2,4,6,8)]

pdf(file="models/biochem-v3.mcmc_posteriors.pdf", width=11, height=8)
par(mar=c(8, 3, 1, 1)+.1)
plot(atvecs, boxes[3,], col=typecols[colvec], pch=20, ylim=range(boxes),
    xaxt='n', xlab='', ylab='parameter value')
axis(1, at=(seq_along(mnames)-1)*(8+xspace) + mean(seq_along(mcmcs)), labels=mnames, las=3)
abline(v=(seq_along(mnames)-1)[-1]*(8+xspace) - xspace/2, lty=3, col='grey')
abline(h=0, lty=3)
segments(x0=atvecs, y0=boxes[1,], y1=boxes[5,], lty=1, col=typecols[colvec])
# segments(x0=atvecs, y0=boxes[2,], y1=boxes[4,], lty=1, col=colvec)
legend("topright", pch=20, col=typecols[seq_along(mcnames)], legend=mcnames)
dev.off()

# zoomed
pdf(file="models/biochem-v3.mcmc_posteriors_zoomed.pdf", width=11, height=8)
plot(atvecs, boxes[3,], col=typecols[colvec], pch=20, ylim=c(-1.2,1.8),
    xaxt='n', xlab='', ylab='parameter value')
axis(1, at=(seq_along(mnames)-1)*(8+xspace) + mean(seq_along(mcmcs)), labels=mnames, las=3)
abline(v=(seq_along(mnames)-1)[-1]*(8+xspace) - xspace/2, lty=3, col='grey')
abline(h=0, lty=3)
segments(x0=atvecs, y0=boxes[1,], y1=boxes[5,], lty=1, col=typecols[colvec])
# segments(x0=atvecs, y0=boxes[2,], y1=boxes[4,], lty=1, col=colvec)
legend("topright", pch=20, col=typecols[seq_along(mcnames)], legend=mcnames)
dev.off()

