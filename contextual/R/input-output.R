###
# input and output functions
#


read.config.counts <- function (infile) {
    count.paramstring <- scan(infile,what='char',nlines=1,sep="\n")
    return( if (substr(count.paramstring,1,1)=="#") {
           jsonlite::fromJSON(gsub("^#*","",count.paramstring),simplifyMatrix=FALSE)
        } else { NULL } )
}

#' Read in a File of Counts 
#' 
#' The file of counts should be of the following form (including the leading commented line):
#'
#'        # { "leftwin" : 3 }
#'        taxon1  taxon2 ... count
#'         AAAA      AA ...    19
#'         CAAA      AA ...     2
#'         GAAA      AA ...     6
#'         TAAA      AA ...     3
#'
#' ... and convert it to a 'tuplecounts' object
#' optionally passing in the orderings of the rows and columns
#'
#' If leftwin is not specified it will be looked for *in json format* on the
#' first line of the file after '#' ... mind the braces and quotes!
#'
#' @return A tuplecounts object.
#'
#' @export
read.counts <- function (infile,leftwin,bases,longpats,shortpats,skip=0) {
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

#' Read in JSON Config File
#'
#' ... or directly from a text string. Automatically fills in zero-length
#' parameters, and fills in selfactors if it isn't there, or moves values from
#' selpats if it is.
#'
#' @return A model configuration (list).
#'
#' @examples
#'
#' model_json <- '
#' {
#'     "bases" : [ "X", "O" ],
#'     "initfreqs" : [ 0.5, 0.5 ],
#'     "mutpats" : [
#'         [ [ "XO", "OX" ] ]
#'     ],
#'     "mutrates" : [ 1 ],
#'     "mutrates.scale" : [ 0.1 ]
#' }
#' '
#' model_config <- read.config(json=model_json)
#'
#' @export
read.config <- function (configfile,quiet=FALSE,json) {
    if (!missing(configfile)&&is.null(configfile)) { cat("Config: NULL.\n"); return(NULL) }
    if (missing(json)) { 
        con <- file(configfile, open="r")
        json <- paste(readLines(con, warn = FALSE), collapse = "\n")
        close(con)
    }
    config <- jsonlite::fromJSON(json,simplifyMatrix=FALSE,simplifyDataFrame=FALSE)
    if (!quiet) { cat("Config: ", jsonlite::toJSON(config), "\n\n") }
    config <- .parse.selpats(config)
    check.config(config)
    return(config)
}

#' Check a configuration has obvious errors.
#'
#' @param config A configuration.
#'
#' @return TRUE or FALSE
check.config <- function (config) {
    ret <- TRUE
    if (!is.null(config$bases) && !is.null(config$initfreqs) && length(config$bases) != length(config$initfreqs)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    if (!is.null(config$mutpats) && !is.null(config$mutrates) && length(config$mutpats) != length(config$mutrates)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    if (!is.null(config$mutpats) && !is.null(config$mutrates.scale) && length(config$mutpats) != length(config$mutrates.scale)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    if (!is.null(config$mutpats) && !is.null(config$mutrates.prior) && length(config$mutpats) != length(config$mutrates.prior)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    if (!is.null(config$selpats) && !is.null(config$selcoef) && length(config$selpats) != length(config$selcoef)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    if (!is.null(config$selpats) && !is.null(config$selcoef.scale) && length(config$selpats) != length(config$selcoef.scale)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    if (!is.null(config$selpats) && !is.null(config$selcoef.prior) && length(config$selpats) != length(config$selcoef.prior)) {
        warning("bases and initfreqs not the same length in configuration")
        ret <- FALSE
    }
    return(ret)
}

#' Fill in Default Values in a Configuration List
#'
#' Fill in default values in a config list based on 'defaults', or length zeros if it makes sense
#'   for mutpats, mutrates, selpats, selfactors, selcoef, bases, fixfn, fixfn.params
#'
#' Note that will NOT fill in default mutrates or selcoef if none are available in defaults
#'   and the corresponding patterns are not length zero
#'   (for instance, models for makegenmat don't need mutrates)
#'
#' @return A model configuration (list).
#'
#' @export
fill.default.config <- function (config, defaults=NULL) {
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

#' Treeify a Configuration.
#'
#' Turn the text representation of a tree (in Newick) into a tree object (as
#' from the ape package), and do associated error checks.  If no tree is
#' specified at all, a stick tree is inserted. Tips and nodes in the tree must be labeled,
#' and branch lengths can be specified in the argument tlen.
#'
#' @return A model configuration (list).
#'
#' @examples
#'
#' # simple two-taxon tree
#' tree_model_json <- '{
#'     "tree" : [ "(sp1 : 0.8, sp2 : 1.2) root;" ],
#'     "bases" : [ "X", "O" ],
#'     "initfreqs" : [ 0.5, 0.5 ],
#'     "sp1" : {
#'         "mutpats" : [ [ [ "XO", "OX" ] ] ],
#'         "mutrates" : [ 1 ]
#'     }, "sp2" : {
#'         "mutpats" : [ [ [ "XO", "OX" ] ] ],
#'         "mutrates" : [ 0.5 ]
#'     } } '
#' orig_config <- read.config(json=tree_model_json)
#' tree_config <- treeify.config(orig_config)
#' # text
#' orig_config$tree
#' # tree object
#' tree_config$tree
#'
#' @export
treeify.config <- function (config,tlen=NULL) {
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

#' Check that A Genmatrix is Compatible with a Configuration
#'
#' @return Returns TRUE if successful; else throws an error.
#'
#' @export
check.genmatrix <- function (config,genmatrix) {
    # TODO: check that fixfn agrees too
    for (xname in c("bases", "mutpats", "selpats", "selfactors")) {
        if (!isTRUE(all.equal( slot(genmatrix,xname), config[[xname]] ))) {
            stop("Precomputed generator matrix does not agree with configuration values for: ", xname)
        }
    }
    return(TRUE)
}

#' Show the names of the nodes of a tree.
#'
#' @export
nodenames <- function (tr) { selfname( c( tr$tip.label, tr$node.label ) ) }

#' @describeIn nodenames Assign node names to a tree.
#' @export
"nodenames<-" <- function (tr,value) { tr$tip.label <- value[seq_along(tr$tip.label)]; tr$node.label <- value[length(tr$tip.label)+seq(1,length.out=tr$Nnode)]; return(tr) }

#' @describeIn nodenames Find the name of a root of a tree.
#' @export
rootname <- function (tr) { nodenames(tr)[ get.root(tr) ] }

#' @describeIn nodenames Assign the name of a root of a tree.
#' @export
"rootname<-" <- function (tr,value) { tr$node.label[get.root(tr)-length(tr$tip.label)] <- value; return(tr) }

#' @describeIn nodenames Return the name of the downstream node of each edge.
#' @export
edge.labels <- function (tr) {
    c(tr$tip.label,tr$node.label)[tr$edge[,2]]
}

#' @describeIn nodenames Return the index of the root in (tips,nodes) order
#' @export
get.root <- function (tr) {
    setdiff( tr$edge[,1], tr$edge[,2] )
}

#' @describeIn nodenames Return index of the parent of node in (tips,nodes) order
#' @export
get.parent <- function (node,tr) {
    tr$edge[match(node,tr$edge[,2]),1]
}

#' @describeIn nodenames Return pairs in node (or NA if none exists)
#' @export
get.cherries <- function (node,tr) {
    parents <- get.parent(node,tr)
    siblings <- outer( parents, parents, "==" )
    sib.indices <- which( siblings & upper.tri(siblings) , arr.ind=TRUE )
    cbind( node[sib.indices[,1]], node[sib.indices[,2]] )
}

####

#' Get a Fixation Function
#'
#' Turn fixfn into an actual function
#'   either by looking it up as a name
#'   or parsing it directly
#' also, check the arguments match fixfn.params.
#'
#' @export
parse.fixfn <- function (fixfn,fixfn.params) {
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

#' Track down references in a config file.
#'
#' Follow pointers in config: if config[[x]] is a string, it refers to another entry.
config.dereference <- function (config, x) {
    sapply(x, function (xx) { n <- 1; while (n < 20 && is.character(config[[xx]]) & (length(config[[xx]])==1) ) { xx <- config[[xx]]; n<-n+1 }; xx } )
}

#' Parse selpats.
#'
#' check if (selpats,selfactors) info is combined into selpats
#'  and separate out if so
#' So: selpats should be EITHER:
#'   - a list of character vectors, OR
#'   - a list of named lists of numeric values
#'
#' @rdname parse.selpats
.parse.selpats <- function (config) {
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

#' Parse models.
#'
#' Check that all models are specified,
#'  fill in defaults,
#'  turn fixfn into functions, etc.
#'
#' @return A model configuration.
#'
#' @examples
#'
#' model_json <- '{
#'     "tree" : [ "(sp1 : 0.8, (sp2 : 1.2, sp3 : 1.0)) root;" ],
#'     "bases" : [ "X", "O" ],
#'     "initfreqs" : [ 0.5, 0.5 ],
#'     "selpats" : [ "X" ],
#'     "selcoef" : [ 0.1 ],
#'     "sp1" : {
#'         "mutpats" : [ [ [ "XO", "OX" ] ] ],
#'         "mutrates" : [ 1 ],
#'         "fixfn" : "popgen.fixfn"
#'     }, "sp2" : {
#'         "mutpats" : [ [ [ "XO", "OX" ] ] ],
#'         "mutrates" : [ 0.5 ]
#'     }, "sp3" : "sp2" } '
#' config <- treeify.config(read.config(json=model_json))
#' config <- parse.models(config)
#' # note that default specified at top level has been copied over
#' config$sp1$selpats
#'
#' @export parse.models
parse.models <- function (config,do.fixfns=TRUE) {
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

# miscellaneous utility functions

#' Name a list or vector with itself.
#' @return Itself, named.
selfname <- function (x) { names(x) <- x; return(x) }
