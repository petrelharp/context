###
# input and output functions
#
# TODO: move some functions over from the other files
###
# frame_files <- lapply(sys.frames(), function(x) x$ofile)
# frame_files <- Filter(Negate(is.null), frame_files)
# .PATH <- dirname(frame_files[[length(frame_files)]])
# source(paste(.PATH,"/helper-fns.R",sep=''))

source("helper-fns.R")

library(jsonlite)

read.config.counts <- function (infile) {
    count.paramstring <- scan(infile,what='char',nlines=1,sep="\n")
    return( if (substr(count.paramstring,1,1)=="#") {
            fromJSON(gsub("^#*","",count.paramstring),simplifyMatrix=FALSE)
        } else { NULL } )
}

read.counts <- function (infile,leftwin,bases,longpats,shortpats,skip=0) {
    # read in a file of counts of the following form:
    #        # { "leftwin" : 3 }
    #        taxon1  taxon2 ... count
    #         AAAA      AA ...    19
    #         CAAA      AA ...     2
    #         GAAA      AA ...     6
    #         TAAA      AA ...     3
    # ... and convert it to a 'tuplecounts' object
    # optionally passing in the orderings of the rows and columns
    #
    # If leftwin is not specified it will be looked for *in json format* on the first line of the file after '#'
    #  ... mind the braces and quotes!
    count.params <- read.config.counts(infile)
    if (!is.null(count.params) && skip==0) { skip <- 1 }
    if (missing(leftwin)) { leftwin <- count.params$leftwin }
    if (is.null(leftwin)) { stop("leftwin not specified.") }
    count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE,skip=skip)
    longwin <- nchar( count.table[1,1] )
    shortwin <- nchar( count.table[1,2] )
    taxa <- colnames(count.table)[-ncol(count.table)]
    if ( missing(bases) ) { bases <- sort( unique( unlist( strsplit( count.table[,1], "" ) ) ) ) }
    if ( missing(longpats) ) { longpats <- getpatterns(longwin,bases) }
    if ( missing(shortpats) ) { shortpats <- getpatterns(shortwin,bases) }
    counts <- Matrix(0,nrow=length(longpats),ncol=length(shortpats)^(length(taxa)-1))
    rownames(counts) <- longpats
    stopifnot( all(count.table[,1] %in% rownames(counts)) )
    stopifnot( all( apply( count.table[,-c(1,ncol(count.table)),drop=FALSE], 2, "%in%", shortpats ) ) )
    colpatterns <- do.call( expand.grid, list( shortpats )[rep.int(1,length(taxa)-1)] )
    colnames(colpatterns) <- taxa[-1]
    colnames(counts) <- apply(colpatterns, 1, paste, collapse='.')
    input.names <- apply( count.table[,-c(1,ncol(count.table)),drop=FALSE], 1, paste, collapse='.' )
    counts[cbind( match(count.table[,1],rownames(counts)), match(input.names,colnames(counts)) )] <- count.table[,ncol(count.table)]
    # stopifnot( sum(counts) == sum(count.table[,ncol(count.table)]) )
    return( new("tuplecounts", 
            counts=counts, 
            leftwin=leftwin, 
            bases=bases,
            colpatterns=colpatterns,
            rowtaxon=taxa[1]
            ) )
}



###
# Config file stuff:

read.config <- function (configfile,quiet=FALSE) {
    # read in JSON config file
    if (is.null(configfile)) { cat("Config: NULL.\n"); return(NULL) }
    con <- openread(configfile)
    json <- paste(readLines(con, warn = FALSE), collapse = "\n")
    close(con)
    config <- fromJSON(json,simplifyMatrix=FALSE,simplifyDataFrame=FALSE)
    if (!quiet) { cat("Config: ", toJSON(config), "\n\n") }
    # fill in zero-length parameters
    # fill in selfactors if it isn't there,
    #  or move values from selpats
    config <- .parse.selpats(config)
    return(config)
}

fill.default.config <- function (config, defaults=NULL) {
    # fill in default values in a config list based on 'defaults', or length zeros if it makes sense
    #   for mutpats, mutrates, selpats, selfactors, selcoef, bases, fixfn, fixfn.params
    #
    # Note that will NOT fill in default mutrates or selcoef if none are available in defaults
    #   and the corresponding patterns are not length zero
    #   (for instance, models for makegenmat don't need mutrates)
    for (x in c("mutpats","selpats","bases","fixfn.params","fixfn.params.scale")) {
        if (is.null(config[[x]])) {
            config[[x]] <- if (is.null(defaults[[x]])) { list() } else { defaults[[x]] }
        }
    }
    if (is.null(config[["selfactors"]])) { 
        config[["selfactors"]] <- if (is.null(defaults[["selfactors"]])) { lapply( config$selpats, sapply, function(x)1.0 ) } else { defaults[["selfactors"]] }
    }
    if (is.null(config[["fixfn"]])) { 
        config[["fixfn"]] <- if(is.null(defaults[["fixfn"]])) { null.fixfn } else { defaults[["fixfn"]] }
    }
    if (length(config$selpats)==0) {
        if (is.null(config$selcoef)) { 
            config$selcoef <- numeric(0) 
        }
        if (is.null(config$selcoef.scale)) {
            config$selcoef.scale <- numeric(0) 
        }
    }
    if (is.null(config$selcoef) && !is.null(defaults$selcoef)) {
        config$selcoef <- defaults$selcoef
    }
    if (is.null(config$mutrates) && !is.null(defaults$mutrates)) {
        config$mutrates <- defaults$mutrates
    }
    if (length(config$selpats)==0 && is.null(config$fixfn)) { 
        config$fixfn <- null.fixfn 
        config$fixfn.params <- list()
    }
    if (length(config$fixfn.params)==0) {
        config$fixfn.params.scale <- list()
    }
    return( config )
}

treeify.config <- function (config,tlen=NULL) {
    # Turn Newick tree representation into a tree, 
    #   and do associated error checks.
    # no tree at all =>  stick tree.
    if (is.null(config$tree)) {
        config <- list( tree="(tip)root;", bases=config$bases, tip=config, initfreqs=config$initfreqs )
    }
    config$tree <- ape::read.tree(text=config$tree)
    if (is.null(config$tree$edge.length)) { 
        # if tree has no edge lengths, bring these in from tlen
        if (!is.null(tlen)) { config$tree$edge.length <- eval(parse(text=tlen)) }
    } else if (!is.null(tlen)) { 
        warning("Branch lengths specified in config file and in tlen; ignoring tlen.")
    }
    if (is.null(config$tree$tip.label) | is.null(config$tree$node.label)) { stop("Please label tips and nodes on the tree.") }
    return(config)
}

check.genmatrix <- function (config,genmatrix) {
    # check that genmatrix is compatible with config
    # TODO: check that fixfn agrees too
    for (xname in c("bases", "mutpats", "selpats", "selfactors")) {
        if (!isTRUE(all.equal( slot(genmatrix,xname), config[[xname]] ))) {
            stop("Precomputed generator matrix does not agree with configuration values for: ", xname)
        }
    }
    return(TRUE)
}

# tree helper functions
nodenames <- function (tr) { selfname( c( tr$tip.label, tr$node.label ) ) }
"nodenames<-" <- function (tr,value) { tr$tip.label <- value[seq_along(tr$tip.label)]; tr$node.label <- value[length(tr$tip.label)+seq(1,length.out=tr$Nnode)]; return(tr) }
rootname <- function (tr) { nodenames(tr)[ get.root(tr) ] }
"rootname<-" <- function (tr,value) { tr$node.label[get.root(tr)-length(tr$tip.label)] <- value; return(tr) }

edge.labels <- function (tr) {
    # return the name of the downstream node of each edge
    c(tr$tip.label,tr$node.label)[tr$edge[,2]]
}

get.root <- function (tr) {
    # return the index of the root in (tips,nodes) order
    setdiff( tr$edge[,1], tr$edge[,2] )
}
get.parent <- function (node,tr) {
    # return index of the parent of node in (tips,nodes) order
    tr$edge[match(node,tr$edge[,2]),1]
}
get.cherries <- function (node,tr) {
    # return pairs in node (or NA if none exists)
    parents <- get.parent(node,tr)
    siblings <- outer( parents, parents, "==" )
    sib.indices <- which( siblings & upper.tri(siblings) , arr.ind=TRUE )
    cbind( node[sib.indices[,1]], node[sib.indices[,2]] )
}

parse.fixfn <- function (fixfn,fixfn.params) {
    # turn fixfn into an actual function
    #   either by looking it up as a name
    #   or parsing it directly
    # also, check the arguments match fixfn.params
    if (is.character(fixfn)) {
        if (exists(fixfn,mode="function")) {
            fixfn <- get(fixfn,mode="function")
        } else {
            fixfn <- eval(parse(text=fixfn))
        }
    }
    if (!missing(fixfn.params)) {
        # check that fixfn.params match what fixfn expects:
        #   first argument is selective differences
        fixfn.argnames <- setdiff(names(as.list(formals(fixfn))),"...")[-1]
        if (!all( fixfn.argnames == names(fixfn.params) )) { 
            stop("fixfn.params (", paste( paste( names(fixfn.params), fixfn.params, sep='=' ), collapse=',' ), ") don't match arguments to fixfn (", paste(fixfn.argnames,collapse=", "), ").")
        }
    }
    return( fixfn )
}

config.dereference <- function (config, x) {
    # follow pointers in config: if config[[x]] is a string, it refers to another entry.
    sapply(x, function (xx) { n <- 1; while (n < 20 && is.character(config[[xx]]) & (length(config[[xx]])==1) ) { xx <- config[[xx]]; n<-n+1 }; xx } )
}

.parse.selpats <- function (config) {
    # check if (selpats,selfactors) info is combined into selpats
    #  and separate out if so
    # So: selpats should be EITHER:
    #   - a list of character vectors, OR
    #   - a list of named lists of numeric values
    if ( !is.null(names(config$selpats)) && is.list(config$selpats) && (!is.list(config$selpats[[1]]) && !is.vector(config$selpats[[1]])) ) { 
        warning("selpats is a named list; did you forget enclosing parentheses?") 
    }
    if (!is.null(config$selpats) && (length(config$selpats)>0)) { 
        config$selpats <- lapply(config$selpats,unlist)
        if (! "selfactors" %in% names(config)) {
            # this takes selpats like c( "XO" = 2.0, "OX" = 1.0 )
            #   and puts the names in selpats and the numbers in selfactors
            #   and makes selfactors a vector of 1.0's otherwise.
            # selfactors retains the names, and so is redundant, but this is not to be relied on.
            config$selfactors <- config$selpats
            no.numbers <- sapply(config$selfactors,is.character)  # these are ones that just get 1.0's
            config$selfactors[no.numbers] <- lapply( config$selfactors[no.numbers], function (x) { y <- rep(1,length(x)); names(y) <- x; return(y) } )
            config$selpats[!no.numbers] <- lapply(config$selpats[!no.numbers], names)
        } else {
            if (any(sapply(config$selpats,is.numeric))) {
                stop("selfactors specified in multiple place.")
            }
        }
    }
    if (!is.null(config$selfactors)) {
        config$selfactors <- lapply( config$selfactors, function (x) { x*1.0 } )  # ensure selfactors is numeric (not integer)
    }
    return(config)
}

parse.models <- function (config,do.fixfns=TRUE) {
    # Check that all models are specified,
    #  fill in defaults,
    #  turn fixfn into functions, etc.
    #
    # Edges are labeled by the node/tip below them:
    nodenames <- nodenames(config$tree)
    rootname <- rootname(config$tree)
    # name root if isn't already (easy to miss)
    if (rootname=="") { rootname <- rootname(config$tree) <- "root"; nodenames <- nodenames(config$tree) }
    edgenodes <- selfname( setdiff( nodenames, rootname ) )
    if (!all(edgenodes %in% names(config)) ) { stop("Must specify named models for each edge in the tree.") }
    # if in config an edge has a string rather than a model, it refers to something else in config
    #  ... so these are the names of the actual models
    edgemodels <- unique( config.dereference(config,edgenodes ) )
    if (any(is.na(edgemodels))) { stop("Must specify named models for each edge in the tree.") }
    config$.models <- edgemodels
    if (!is.null(config$.models)) { names(config$.models) <- edgemodels }
    # parse selpats into selpats,selfactors, if these are there
    config <- .parse.selpats(config)
    # get the fixfns in there
    for (modname in edgemodels) {
        # parse selpats into selpats,selfactors
        config[[modname]] <- .parse.selpats(config[[modname]])
        # put in defaults e.g. for no selection
        config[[modname]] <- fill.default.config(config[[modname]],defaults=config)
        if (do.fixfns) {  # should come last
            # turn fixfn into a function and check we have the right parameters
            config[[modname]]$fixfn <- parse.fixfn(config[[modname]]$fixfn,config[[modname]]$fixfn.params)
        }
        stopifnot( with( config[[modname]], ( length(mutpats) == length(mutrates) ) && ( is.null(config[[modname]]$selcoef) || ( ( length(selpats) == length(selcoef) ) && ( length(selcoef) == length(selfactors) ) ) ) ) )
    }
    return(config)
}


## plotting functions

likelihood.surface <- function (model, tlen=1, 
        plot.it=TRUE, ask=FALSE, progress=TRUE,
        mutrates=model@mutrates, selcoef=model@selcoef, 
        do.parallel=("parallel" %in% .packages(all.available=TRUE)),
        numcores=if (do.parallel) { parallel::detectCores() } else { 1 },
        ngrid=10) {
    # Compute the likelihood on a grid of points about the parameters in 'model'
    # and plot these (note: plots more than one figure)
    this.lapply <- if ( do.parallel && numcores>1 ) { function (...) { parallel::mclapply( ..., mc.cores=numcores ) } } else { lapply }

    timevec <- c( rep(as.numeric(tlen),nmuts(model)), rep(1,length(coef(model))-nmuts(model)) )
    if (!missing(mutrates)) { model@mutrates <- mutrates*tlen }
    if (!missing(selcoef)) { model@selcoef <- selcoef }
    pvec <- coef(model)/timevec
    f <- function (x) { model@likfun((x*timevec)[model@results$use.par]) }

    # contours
    jk <- combn(which(model@results$use.par),2)
    grids <- this.lapply( 1:ncol(jk), function (ii) {
            j <- jk[1,ii]
            k <- jk[2,ii]
            inds <- 1:length(pvec)
            inds[j] <- 1+length(pvec)
            inds[k] <- 2+length(pvec)
            pvals <- seq( .9*min(pvec[j]), 1.1*max(pvec[j]), length.out=ngrid )
            qvals <- seq( .9*min(pvec[k]), 1.1*max(pvec[k]), length.out=ngrid )
            lvals <- matrix( NA, nrow=length(pvals), ncol=length(qvals) )
            for (jj in seq_along(pvals)) {
                for (kk in seq_along(qvals)) {
                    lvals[jj,kk] <- f( c(pvec,pvals[jj],qvals[kk])[inds] )
                }
            }
            if (progress) { cat(".") }
            glist <- list( pvals, qvals, lvals )
            names( glist ) <- c( names(coef(model)[model@results$use.par])[c(j,k)], "lvals" )
            return(glist)
        } )
    if (plot.it) {
        for (ii in seq_along(grids)) {
            j <- jk[1,ii]
            k <- jk[2,ii]
            cols <- colorspace::diverge_hcl(64)
            plot( grids[[ii]][[1]][row(grids[[ii]]$lvals)], grids[[ii]][[2]][col(grids[[ii]]$lvals)], cex=3, pch=20, 
                    col=cols[cut(grids[[ii]]$lvals,breaks=length(cols))], ask=ask,
                    xlab=names(coef(model)[model@results$use.par])[j],
                    ylab=names(coef(model)[model@results$use.par])[k]
                )
            contour( grids[[ii]][[1]], grids[[ii]][[2]], grids[[ii]]$lvals, add=TRUE )
        }
    }
    return( grids )
}
