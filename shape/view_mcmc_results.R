library(contextual)

for (modelname in c("biochem-v5", "biochem-v4", "biochem-v3")) {
    mcmcs <- list.files("RegulatoryFeature-regions-from-axt/", 
                paste0(modelname,"-fit-mcmc.*RData"), recursive=TRUE, full.names=TRUE)

    source("models/mutrate_names.R") # provides mutrate.names

    mcmcs <- mcmcs[sapply(mcmcs, function (mcfile) { load(mcfile); (dim(model@results[["batch"]])[1]>1000) } )]

    # check all the same
    all.mnames <- lapply(mcmcs, function (mcfile) { load(mcfile);
                         gsub(".*[.]","",names(coef(model))[model@results$use.par]) } )
    for (k in seq_along(all.mnames)[-1]) { stopifnot(all(all.mnames[[k]]==all.mnames[[1]])) }

    dir.create(file.path("models",paste0(modelname,"_mcmc_traces")), recursive=TRUE, showWarnings=FALSE)
    for (mcfile in mcmcs) {
        load(mcfile)
        mnames <- names(mutrate.names)[match(gsub(".*[.]","",names(coef(model))[model@results$use.par]), 
                                             mutrate.names)]
        png(file=file.path("models", paste0(modelname,"_mcmc_traces"), 
                           paste0(gsub("/", "_", gsub("R[^/]*/*","",gsub(".RData","",mcfile))), "_trace.png")),
            width=11*300, height=8*300, res=300, type='cairo')
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
        dev.off()
    }

    for (nonnq in c(TRUE, FALSE)) {
        nameq <- if (nonnq) { "nonnegative" } else { "unconstrained" }

        these_mcmcs <- mcmcs[sapply(mcmcs, function (mcfile) { load(mcfile); (dim(model@results[["batch"]])[1]>20000) && (all(model@results[["batch"]]>=0) == nonnq) } )]
        if (length(these_mcmcs) > 0) {

            xspace <- 4
            boxes <- do.call(cbind, lapply(lapply(these_mcmcs, function (mcfile) { load(mcfile); boxplot(as.data.frame(model@results[["batch"]][seq(20001,40000),]), plot=FALSE) } ), "[[", "stats"))
            write.table(boxes, file=sprintf("models/%s_%s_mcmc_posteriors.tsv", modelname,nameq), row.names=FALSE)

            atvecs <- unlist(lapply(seq_along(these_mcmcs), function (k) { (seq_along(mnames)-1)*(8+xspace) + k }))
            colvec <- rep(1:8, each=length(mnames))
            typecols <- RColorBrewer::brewer.pal("Paired", n=8)[c(1,3,5,7,2,4,6,8)]
            mcnames <- gsub("noOverlap", "nongene", gsub("overlap", "gene", gsub("-knownGeneTx", "", sapply(these_mcmcs, function (mcfile) gsub("/", " : ", gsub("R[^/]*[/]*","",dirname(dirname(mcfile))) )))))


            ylims <- list(range(boxes), "_zoomed"=c(-1.2,1.8))
            for (k in seq_along(ylims)) {
                pdf(file=sprintf("models/%s_%s_mcmc_posteriors%s.pdf",modelname,nameq,names(ylims)[k]), width=6.5, height=4, pointsize=10)
                layout(t(1:2), widths=c(3,1))
                par(mar=c(6, 4, 1, 1)+.1)
                plot(atvecs, boxes[3,], col=typecols[colvec], pch=20, ylim=ylims[[k]],
                    xaxt='n', xlab='', ylab='parameter value')
                axis(1, at=(seq_along(mnames)-1)*(8+xspace) + mean(seq_along(these_mcmcs)), labels=mnames, las=3)
                abline(v=(seq_along(mnames)-1)[-1]*(8+xspace) - xspace/2, lty=3, col='grey')
                abline(h=0, lty=3)
                segments(x0=atvecs, y0=boxes[1,], y1=boxes[5,], lty=1, col=typecols[colvec])
                # segments(x0=atvecs, y0=boxes[2,], y1=boxes[4,], lty=1, col=colvec)
                par(mar=c(0, 0, 1, 1)+.1)
                plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
                legend("topright", pch=20, col=typecols[seq_along(mcnames)], legend=mcnames, cex=0.5, bty='n')
                dev.off()
            }
        }
    }
}
