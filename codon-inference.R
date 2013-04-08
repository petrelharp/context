#!/usr/bin/R
source("codon-inference-fns.R")  # list of codons, etc 

# size of window on either side of the focal site
lwin <- rwin <- 2
winlen <- lwin+rwin+1

patterns <- do.call( expand.grid, rep( list(bases), winlen ) )
npatt <- nrow(patterns)
baseind <- nbases^(0:(winlen-1))
rownames(patterns) <- apply(patterns,1,paste,collapse="")


# list of patterns that have mutation rates
mutpats <- c(
        combn( bases, 2, simplify=FALSE ),  # single-base rates
        list( c("CG","TG"), c("CG","CA") ), # CpG rates
        # do.call( c, with( subset(codons, aa%in%synons), tapply( levels(codon)[as.numeric(codon)], droplevels(aa), function(x)combn(x, 2, simplify=FALSE) ) ) ),  # synonymous codon rates 
    NULL )


# list of patterns with selection coefficients
selpats <- c(
        codons$codon[codons$aa %in% synons],
    NULL )

# 
