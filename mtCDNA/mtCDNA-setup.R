#!/usr/bin/Rscript --vanilla
# setup for mtCDNA inference

source("../codons.R")
codonpairs <- combn( as.character(codons$codon), 2 )
ndiffs <- apply( codonpairs, 2, function (x) { sum( do.call( "!=",  strsplit(x,"") ) ) } )
codonpairs <- codonpairs[,ndiffs==1]

mutpats <- c(
        apply(codonpairs,2,list),  
        apply(codonpairs[2:1,],2,list), 
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate, why not?
    ) 
# selpats <- as.list( as.character(codons$codon) )
# fixfn <- popgen.fixfn  # takes Ne as additional parameter
selpats <- list()
fixfn <- function (...) { 1 }

mutrates <- rep(1e-8,length(mutpats))
selcof <- numeric(0)
