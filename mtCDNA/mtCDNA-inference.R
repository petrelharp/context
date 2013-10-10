#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tree-cpg.R .\
"

option_list <- list(
        make_option( c("-c","--infile"), type="character", default="mtCDNApri-sub.nuc", help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-v","--tprior"), type="double", default=.5, help="Parameter for Beta prior on branch length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-e","--countfile"), type="character", default="", help="Record tuple counts in this file. [default: countdata/counts-(lwin)-(win)-(rwin).RData]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

winlen <- lwin+win+rwin

if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') }

if (is.null(infile)) { cat("Run\n  mtCDNA-inference.R -h\n for help.") }

if (logfile!="") {
    logfile <- gsub(".RData",".Rout",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message", split=interactive()) 
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}

if (countfile=="") { countfile <- paste(paste("countdata/counts",lwin,win,rwin,sep="-"),".RData",sep="") }


# took mtCDNApri.nuc from paml/examples/mtCDNA
#  and: removed extra spaces on first line (only one space between numbers)
#       made each sequence be on a single line
#       removed trailing text

require(Biostrings)

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

basedir <- gsub(".nuc","",infile,fixed=TRUE)
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')


if (file.exists(gmfile)) {
    load(gmfile)
} else {
    if (meanboundary>0) {
        genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen[1], selcoef=numeric(0), boundary=boundary )
    } else {
        genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen[1], selcoef=numeric(0), boundary=boundary )
    }
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

# all pairwise counts
if (!file.exists(countfile)) {
    mtCDNA <- as( readDNAMultipleAlignment(infile, format="phylip" ), "DNAStringSet" )
    mtCDNA <- mtCDNA[ c("human", "chimpanzee", "bonobo", "gorilla") ]

    require(parallel)
    counts <- mclapply( names( mtCDNA ), function (x) {
            x <- lapply( names( mtCDNA ), function (y) {
                    counttrans( rownames(projmatrix), colnames(projmatrix), mtCDNA[[x]],  mtCDNA[[y]],  lwin=lwin, shift=3 ) 
                } ) 
            names(x) <- names( mtCDNA )
            return(x)
        }, mc.cores=6 )
    names(counts) <- names( mtCDNA )
    save( counts, file=countfile )
} else {
    load( countfile )
}

# not used in likelihood:
initcounts <- lapply( lapply(counts,"[[",1), lapply, rowSums )  # initial counts for (x,y) only depends on x

# set up root distribution
nmuts <- length(mutpats)
nfreqs <- length(bases)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
patcomp <- t(patcomp) # column-wise more efficient

which.taxa <- c("human","chimpanzee")
# which.frame will denote which element of counts[[x]][[y]] to use

# composite likelihood with taxa x->y + y->x
likfun <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:nfreqs)]
    initfreqs <- initfreqs/sum(initfreqs)
    patfreqs <- log(initfreqs)[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- exp( colSums( patfreqs ) )
    # these are collapsed transition matrix
    updownbranch <- list(  # note "up" branch is from simpler summaries
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) ),
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        )
    # if (any(sapply(updownbranch,function(x) any(!is.numeric(x))|any(x<0)))) { browser() }
    # return negative log-likelihood plus a penalty to keep initfreqs summing to (almost) 1
    return( 
                (-1) * ( sum( counts[[which.taxa[[1]]]][[which.taxa[[2]]]][[which.frame]] * log(updownbranch[[1]]) ) + sum( counts[[which.taxa[[2]]]][[which.taxa[[1]]]][[which.frame]] * log(updownbranch[[2]]) ) ) 
                + 100*(sum(initfreqs)-1)^2
            )
}

initparams <- c( rel.tlen=0.5, 6e6*rep(1e-8,nmuts)/30, rep(1/nfreqs,nfreqs) )  # reasonable for hu-ch?
stopifnot( is.finite( likfun(initparams) ) )

save( lwin, rwin, win, winlen, boundary, meanboundary, gmfile, projmatrix, subtransmatrix, counts, initcounts, file=countfile )

