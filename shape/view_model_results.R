library(contextual)

all.results <- read.table("model-selection-results.tsv", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
# remove non-strand-symmetric model
all.results <- subset(all.results, model!="base-model")
all.results <- all.results[,!apply(is.na(all.results),2,all)]
# fix up wierd duplicated columns
na.sum <- function (x,y) { ifelse(is.na(x), ifelse(is.na(y), NA, y), ifelse(is.na(y), x, (x+y)/2)) }
all.results[["fit:reference.CG->CA|CG->TG"]] <- na.sum(all.results[["fit:reference.CG->CA|CG->TG"]], all.results[["fit:reference.CG->TG|CG->CA"]])
all.results[["fit:derived.CG->CA|CG->TG"]] <- na.sum(all.results[["fit:derived.CG->CA|CG->TG"]], all.results[["fit:derived.CG->TG|CG->CA"]])
all.results[["fit:reference.CG->CT|CG->AG"]] <- na.sum(all.results[["fit:reference.CG->CT|CG->AG"]], all.results[["fit:reference.CG->CT|CG->AG.1"]])
all.results <- all.results[setdiff(names(all.results), c("fit:reference.CG->TG|CG->CA","fit:derived.CG->TG|CG->CA","fit:reference.CG->CT|CG->AG.1"))]



all.info <- all.results[!grepl("^fit:", names(all.results))][,-(1:2)]
all.info$type <- factor(all.info$type)
all.info$model <- factor(all.info$model)
# color will be regulatory type
all.info$color <- as.numeric(all.info$type)
all.info$cex <- all.results$longwin/4
all.info$pch <- 1 + all.info$overlap
all.info$model.color <- rainbow(nlevels(all.info$model))[as.numeric(all.info$model)]

source("models/mutrate_names.R")

mutrate.cols <- grepl("->", names(all.results))
mutrates.branch <- all.results[,mutrate.cols]

mutrates.info <- data.frame( param=gsub("fit.","",names(all.results)[mutrate.cols]), stringsAsFactors=FALSE, check.names=FALSE )
mutrates.info$reference <- grepl("^reference", mutrates.info$param)
mutrates.info$mutpat <- gsub("[^.]*\\.", "", mutrates.info$param)
mutrates.info$name <- factor(names(mutrate.names)[match(mutrates.info$mutpat, mutrate.names)], levels=names(mutrate.names))

# lapply(mutrate.names, function (x) { match(paste0("fit:", c("reference", "derived"), ".", x), names(all.results)) } )

mutrates <- do.call(cbind, lapply(mutrate.names, function (x) {
        mcols <- intersect(paste0("fit:", c("reference", "derived"), ".", x), names(all.results))
        if (length(mcols)>0) {
            out <- all.results[,mcols[1]]
            if (length(mcols)>1) {
                out <- na.sum(out, all.results[,mcols[2]])
            }
        } else { out <- NULL }
        return(out)
    } ) )

# call out those without single-base rates
singlerates <- unique(unlist(list(c("A<->T","C<->G", "A->C", "A->G", "C->A", "C->T"), c("A<->T", "A<->C", "A<->G", "C<->G"))))
singlewierd <- (rowSums(mutrates[,singlerates]==1e-6, na.rm=TRUE) > 0)


plot_fun <- function (colorvec, colors, labels, subset=TRUE) {
    # color by regulatory region type
    layout(1:2, heights=c(1,1.3))
    par(mar=c(0, 3, 1, 1)+.1)
    plot(jitter(as.vector(col(mutrates[subset,]))),
         as.vector(mutrates[subset,]),
         col=colorvec[subset], 
         cex=all.info$cex[subset],
         pch=all.info$pch[subset],
         xaxt='n', xlab='', ylab='mutation rate')
    legend("topright", col=colors, pch=1, legend=labels)
    par(mar=c(10, 3, 0, 1)+.1)
    plot(jitter(as.vector(col(mutrates[subset,]))),
         as.vector(mutrates[subset,]), log='y',
         col=colorvec[subset], 
         cex=all.info$cex[subset],
         pch=all.info$pch[subset],
         xaxt='n', xlab='', ylab='mutation rate')
    axis(1, at=1:ncol(mutrates), labels=colnames(mutrates), las=2)
    legend("bottomright", pch=1:2, legend=c("overlap", "noOverlap"))
}

pdf(file="all-model-results.pdf", width=11, height=8, pointsize=10)

plot_fun(all.info$color, 1:nlevels(all.info$type), levels(all.info$type), subset=!singlewierd)

plot_fun(all.info$color, 1:nlevels(all.info$type), levels(all.info$type))

plot_fun(all.info$model.color, rainbow(nlevels(droplevels(all.info$model[!singlewierd]))), 
         levels(droplevels(all.info$model[!singlewierd])), subset=!singlewierd)

plot_fun(all.info$model.color, rainbow(nlevels(all.info$model)), levels(all.info$model))

dev.off()
