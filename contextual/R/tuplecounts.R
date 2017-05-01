
#' tuplecounts Class
#'
#' tuplecounts is a (2D, flattened) contingency table
#'
#' @name tuplecounts-class
#' @rdname tuplecounts-class
#' @exportClass tuplecounts
setClass("tuplecounts",representation(
            leftwin="numeric",
            counts="Matrix",
            bases="character",
            rowtaxon="character",
            colpatterns="data.frame"),
        prototype=list(rowtaxon="long")
    )

# extractor functions
#' @name rowtaxon
#' @rdname tuplecounts-methods
#' @exportMethod rowtaxon
setGeneric("rowtaxon", function(x) { standardGeneric("rowtaxon") })
#' @rdname tuplecounts-methods
#' @aliases rowtaxon,tuplecounts-method
setMethod("rowtaxon", signature=c(x="tuplecounts"), definition=function (x) { x@rowtaxon } )
#' @name coltaxa
#' @rdname tuplecounts-methods
#' @exportMethod coltaxa
setGeneric("coltaxa", function(x) { standardGeneric("coltaxa") })
#' @rdname tuplecounts-methods
#' @aliases coltaxa,tuplecounts-method
setMethod("coltaxa", signature=c(x="tuplecounts"), definition=function (x) { colnames(x@colpatterns) } )
#' @name colpatterns
#' @rdname tuplecounts-methods
#' @exportMethod colpatterns
setGeneric("colpatterns", function(x) { standardGeneric("colpatterns") })
#' @rdname tuplecounts-methods
#' @aliases colpatterns,tuplecounts-method
setMethod("colpatterns", signature=c(x="tuplecounts"), definition=function (x) { x@colpatterns } )
#' @name counts
#' @rdname tuplecounts-methods
#' @exportMethod counts
setGeneric("counts", function(x) { standardGeneric("counts") })
#' @rdname tuplecounts-methods
#' @aliases counts,tuplecounts-method
setMethod("counts", signature=c(x="tuplecounts"), definition=function (x) { x@counts } )
#'
#' This method repackages the counts in a more friendly-looking data frame
#'
#' @name countframe
#' @rdname tuplecounts-methods
#' @exportMethod countframe
setGeneric("countframe", function(x) { standardGeneric("countframe") })
#' @rdname tuplecounts-methods
#' @aliases countframe,tuplecounts-method
setMethod("countframe", signature=c(x="tuplecounts"), definition=function (x) {
        cf <- cbind(
                data.frame( rep.int(rownames(x@counts),ncol(x@counts)) ),
                colpatterns(x)[ rep(1:nrow(colpatterns(x)),each=nrow(x@counts)), ],
                as.numeric(x@counts)
            )
        colnames(cf) <- c( rowtaxon(x), coltaxa(x), "count" )
        return(cf)
    } )

setGeneric("dim")
setGeneric("dimnames")
# things to make tuplecounts act like the matrix inside of it:
#' @rdname tuplecounts-methods
#' @aliases dim,tuplecounts-method
setMethod("dim", signature=c(x="tuplecounts"), definition=function (x) { dim(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases dimnames,tuplecounts-method
setMethod("dimnames", signature=c(x="tuplecounts"), definition=function (x) { dimnames(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases dimnames<-,tuplecounts-method
setMethod("dimnames<-", signature=c(x="tuplecounts",value="ANY"), definition=function (x,value) { dimnames(x@counts)<-value } )
#' @rdname tuplecounts-methods
#' @aliases as.matrix,tuplecounts-method
setMethod("as.matrix", signature=c(x="tuplecounts"), definition=function (x) { as.matrix(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases as.vector,tuplecounts-method
setMethod("as.vector", signature=c(x="tuplecounts"), definition=function (x) { as.vector(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases head,tuplecounts-method
setMethod("head", signature=c(x="tuplecounts"), definition=function (x) { head(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases image,tuplecounts-method
setMethod("image", signature=c(x="tuplecounts"), definition=function (x) { image(x@counts) } )
setMethod("%*%", signature=c(x="tuplecounts",y="ANY"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="ANY",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
# add more fluff to avoid "Note: method with signature ‘ANY#tuplecounts’ chosen for function ..." messages
setMethod("%*%", signature=c(x="tuplecounts",y="tuplecounts"), definition=function (x,y) { x@counts %*% y@counts } )
setMethod("%*%", signature=c(x="dgCMatrix",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("%*%", signature=c(x="tuplecounts",y="dgCMatrix"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="dgTMatrix",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("%*%", signature=c(x="tuplecounts",y="dgTMatrix"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="TsparseMatrix",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("%*%", signature=c(x="tuplecounts",y="TsparseMatrix"), definition=function (x,y) { x@counts %*% y } )
