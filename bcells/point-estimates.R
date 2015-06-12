#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from paired counts file.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", default=NULL, help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of left-hand context." ),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-m","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
print(opt) # this will go in the pbs log
attach(opt)
options(error = quote({dump.frames(to.file = TRUE); q()}))


if (!file.exists(infile)) { stop("Cannot read file ", infile) }

gmname <- gsub(".RData","",basename(gmfile))
basename <- paste(basedir,"/", basename(infile), "-", gmname, "-", jobid, sep='')
datafile <- paste( basename ,"-results.RData",sep='')

logfile <- paste(basename,".Rout",sep='')
logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
sink(file=logcon, type="message")
sink(file=logcon, type="output", split=interactive())   # send both to log file

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''),chdir=TRUE)

date()
cat("basename: ", basename, "\n")
print(opt)


if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop("Can't find generator matrix in ", gmfile, " -- provide file name exactly?")
}

# read in counts (produced with count-paired-tuples.py)
count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)

# check window sizes match
longwin <- nchar(rownames(genmatrix)[1])
stopifnot( longwin == nchar( count.table$reference[1] ) )
rightwin <- longwin - leftwin - nchar( count.table$derived[1] )
shortwin <- longwin - leftwin - rightwin
stopifnot( shortwin == nchar( count.table$derived[1] ) & rightwin > 0 & shortwin > 0 & leftwin > 0 )

# projection matrix
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )

# parse counts into a matrix
counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
rownames(counts) <- rownames(genmatrix)
colnames(counts) <- colnames(projmatrix)
stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
initcounts <- rowSums(counts)

# ad-hoc estimate
adhoc <- countmuts(counts=counts,mutpats=mutpats,leftwin=leftwin)
adhoc <- adhoc[1,]/adhoc[2,]
stopifnot( all( is.finite( adhoc ) )

nmuts <- length(genmatrix@mutpats)
nsel <- length(genmatrix@selpats)
stopifnot(nsel==0)
# (quasi)-likelihood function using all counts -- multinomial
likfun <- function (params) {
    # params are: mutrates*tlen, shape
    genmatrix@x <- update(genmatrix,mutrates=params[1:nmuts],selcoef=numeric(0))
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return NEGATIVE log-likelihood 
    ans <- (-1) * sum( counts * log(subtransmatrix) ) # + (params[length(params)]-1)^2
    if (!is.finite(ans)) { print(params) }
    return(ans)
}


# point estimates
initpar <- adhoc # ,1 )
lbs <- rep(1e-6,nmuts) #, .01 )
ubs <- rep(.5,nmuts) #, 5 )
parscale <- 1e-3 * rep(mean(adhoc),length(adhoc)) #, 1 )

baseval <- likfun(initpar)
stopifnot( is.finite(baseval) )
optim.point.estimate <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(fnscale=abs(baseval), parscale=parscale, maxit=100) )
if (optim.point.estimate$convergence!=0) { warning("optim failed to converge.") }
point.estimate <- optim.point.estimate$par
names(point.estimate) <- mutnames(genmatrix@mutpats)

date()
cat("done with computation.\n")
cat("saving to: ", datafile, "\n")

save( opt, gmname, infile, adhoc, likfun, nmuts, counts, initpar, lbs, ubs, parscale, optim.point.estimate, point.estimate,
     longwin, rightwin, leftwin, shortwin, 
     file=datafile )

print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
sink(NULL); flush(logcon)
