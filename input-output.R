###
# input and output functions
#
# TODO: move some functions over from the other files
###
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
.PATH <- dirname(frame_files[[length(frame_files)]])
source(paste(.PATH,"/helper-fns.R",sep=''))

require(jsonlite)

read.config.counts <- function (infile) {
    count.paramstring <- scan(infile,what='char',nlines=1,sep="\n")
    return( if (substr(count.paramstring,1,1)=="#") {
            fromJSON(gsub("^#*","",count.paramstring),simplifyMatrix=FALSE)
        } else { NULL } )
}

read.counts <- function (infile,leftwin,bases,longpats,shortpats,skip=0) {
    # read in a file of counts of the following form:
    #     # leftwin=1
    #        taxon1  taxon2 ... count
    #         AAAA      AA ...    19
    #         CAAA      AA ...     2
    #         GAAA      AA ...     6
    #         TAAA      AA ...     3
    # ... and convert it to a 'tuplecounts' object
    # optionally passing in the orderings of the rows and columns
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
    config <- fromJSON(json,simplifyMatrix=FALSE)
    if (!quiet) { cat("Config: ", toJSON(config), "\n\n") }
    # fill in zero-length parameters
    if (is.null(config$selpats)) { config$selpats <- list() }
    if (length(config$selpats)==0 && is.null(config$selcoef)) { 
        config$selcoef <- numeric(0) 
        config$selcoef.scale <- numeric(0) 
    }
    if (length(config$selpats)==0 && is.null(config$fixfn)) { 
        config$fixfn <- null.fixfn 
        config$fixfn.params <- list()
    }
    if (length(config$fixfn.params)==0) {
        config$fixfn.params.scale <- list()
    }
    return(config)
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

fill.default.config <- function (config, defaults=NULL) {
    # fill in default values in a config list
    # for selpats, bases, fixfn, fixfn.params
    for (x in c("selpats","bases","fixfn.params")) {
        if (is.null(config[[x]])) {
            config[[x]] <- if (is.null(defaults[[x]])) { list() } else { defaults[[x]] }
        }
    }
    if (is.null(config[["fixfn"]])) { 
        config[["fixfn"]] <- if(is.null(defaults[["fixfn"]])) { null.fixfn } else { defaults[["fixfn"]] }
    }
    return( config )
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
    if (is.null(fixfn)) { fixfn <- null.fixfn }
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
            stop("fixfn.params (", paste( paste( names(fixfn.params), fixfn.params, sep='=' ), collapse=',' ), ") don't match arguments to fixfn (", fixfn.argnames, ").")
        }
    }
    return( fixfn )
}

config.dereference <- function (config, x) {
    # follow pointers in config: if config[[x]] is a string, it refers to another entry.
    sapply(x, function (xx) { n <- 1; while (n < 20 && is.character(config[[xx]]) & (length(config[[xx]])==1) ) { xx <- config[[xx]]; n<-n+1 }; xx } )
}

parse.models <- function (config,do.fixfns=TRUE) {
    # Check that all models are specified, and turn each fixfn into a function.
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
    # get the fixfns in there
    if (do.fixfns) for (modname in edgemodels) {
        # turn fixfn into a function and check we have the right parameters
        config[[modname]]$fixfn <- parse.fixfn(config[[modname]]$fixfn,config[[modname]]$fixfn.params)
        # put in defaults if no selection
        if (is.null(config[[modname]]$selpats)) { config[[modname]]$selpats <- list(); config[[modname]]$selcoef <- numeric(0); config[[modname]]$fixfn <- null.fixfn; config[[modname]]$fixfn.params=list() }
        stopifnot( with( config[[modname]], ( length(mutpats) == length(mutrates) ) && ( length(selpats) == length(selcoef) )))
    }
    return(config)
}
