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

mutrate.names <- c(
        "A<->T" = "A->T|T->A",
        "A<->C" = "A->C|T->G|C->A|G->T",
        "A<->G" = "A->G|T->C|C->T|G->A",
        "C<->G" = "C->G|G->C",
        "A->C" = "A->C|T->G",
        "A->G" = "A->G|T->C",
        "C->A" = "C->A|G->T",
        "C->T" = "C->T|G->A",
        "CpG" = "CG->CA|CG->TG",
        "CpG+UV" = "TCG->TTG|CCG->CTG|CGA->CAA|CGG->CAG",
        "AID short" = "AC->AT|GC->GT|GT->AT|GC->AC",
        "AID WRCY" = "AACT->AATT|AACC->AATC|AGCT->AGTT|AGCC->AGTC|TACT->TATT|TACC->TATC|TGCT->TGTT|TGCC->TGTC|AGTT->AATT|GGTT->GATT|AGCT->AACT|GGCT->GACT|AGTA->AATA|GGTA->GATA|AGCA->AACA|GGCA->GACA",
        "APOBEC" = "TCA->TTA|TCT->TTT|TCA->TGA|TCT->TGT|TGA->TAA|AGA->AAA|TGA->TCA|AGA->ACA",
        "iota double" = "AA->GT|TT->AC",
        "iota single" = "AA->AG|CT->TT",
        "eta" = "TA->TG|TA->CA",
        "CG->CT" = "CG->CT|CG->AG",
        "CG->TG" = "CG->TG|CG->CA",
        "AA->TT" = "AA->TT|TT->AA",
        "GT->AT" = "GT->AT|AC->AT",
        "TC->TT" = "TC->TT|GA->AA",
        "CC->CT" = "CC->CT|GG->AG",
        "CC->TT" = "CC->TT|GG->AA",
        "TT->TX" = "TT->TA|AA->TA|TT->TC|AA->GA|TT->TG|AA->CA|TT->AT|AA->AT|TT->CT|AA->AG|TT->GT|AA->AC",
        "CC->CX" = "CC->CA|GG->TG|CC->CT|GG->AG|CC->CG|GG->CG|CC->AC|GG->GT|CC->TC|GG->GA|CC->GC|GG->GC",
        "CT->CX" = "CT->CA|AG->TG|CT->CC|AG->GG|CT->CG|AG->CG|TC->AC|GA->GT|TC->CC|GA->GG|TC->GC|GA->GC",
        "TT->XX" = "TT->AA|AA->TT|TT->AC|AA->GT|TT->AG|AA->CT|TT->CA|AA->TG|TT->CC|AA->GG|TT->CG|AA->CG|TT->GA|AA->TC|TT->GC|AA->GC|TT->GG|AA->CC",
        "CC->XX" = "CC->AA|GG->TT|CC->AT|GG->AT|CC->AG|GG->CT|CC->TA|GG->TA|CC->TT|GG->AA|CC->TG|GG->CA|CC->GA|GG->TC|CC->GT|GG->AC|CC->GG|GG->CC",
        "CT->XX" = "CT->AA|AG->TT|CT->AC|AG->GT|CT->AG|AG->CT|CT->TA|AG->TA|CT->TC|AG->GA|CT->TG|AG->CA|CT->GA|AG->TC|CT->GC|AG->GC|CT->GG|AG->CC",
        "TC->XX" = "TC->AA|GA->TT|TC->AT|GA->AT|TC->AG|GA->CT|TC->CA|GA->TG|TC->CT|GA->AG|TC->CG|GA->CG|TC->GA|GA->TC|TC->GT|GA->AC|TC->GG|GA->CC",
        "CT->CC" = "CG->CC|CG->GG",
        "TCC->TTC" = "TCC->TTC|AGG->AAG",
        "AT->GT" = "AT->GT",
        "GG->AG" = "GG->AG"
        )

mutrate.cols <- grepl("->", names(all.results))
mutrates.branch <- all.results[,mutrate.cols]

mutrates.info <- data.frame( param=gsub("fit.","",names(all.results)[mutrate.cols]), stringsAsFactors=FALSE, check.names=FALSE )
mutrates.info$reference <- grepl("^reference", mutrates.info$param)
mutrates.info$mutpat <- gsub("[^.]*\\.", "", mutrates.info$param)
mutrates.info$name <- factor(names(mutrate.names)[match(mutrates.info$mutpat, mutrate.names)], levels=names(mutrate.names))

lapply(mutrate.names, function (x) { match(paste0("fit:", c("reference", "derived"), ".", x), names(all.results)) } )

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
singlewierd <- (rowSums(mutrates[,z]==1e-6, na.rm=TRUE) > 0)


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
