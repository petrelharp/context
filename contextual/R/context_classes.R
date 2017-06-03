#' @include genmatrix.R tuplecounts.R
NULL

#' A model of mutation; one for each branch.
#'
#' @name contextModel-class
#' @rdname contextModel-class
#' @exportClass contextModel
setClass("contextModel",
         representation(
                       genmatrix="genmatrix",
                       projmatrix="Matrix",
                       mutrates="numeric",
                       selcoef="numeric",
                       params="numeric"
                    )
            )

#' A model of mutation
#'
#' @name context-class
#' @rdname context-class
#' @exportClass context
setClass("context",
         representation(
                       counts="tuplecounts",
                       likfun="function",
                       results="list",
                       invocation="character"
                    ),
            contains="contextModel"
         )

#' The results of a contextual MCMC
#'
#' @name contextMCMC-class
#' @rdname contextMCMC-class
#' @exportClass contextMCMC
setClass("contextMCMC", 
        representation(
                       mutrates.prior="numeric",
                       selcoef.prior="numeric",
                       fixfn.params.prior="numeric"
                   ),
             contains="context"
         )

setOldClass("phylo")  # phylo in ape is an S3 class

#' A context model on a tree
#'
#' @name contextTree-class
#' @rdname contextTree-class
#' @exportClass contextTree
setClass("contextTree", representation(
                       counts="tuplecounts",
                       tree="phylo",  # as in 'ape'
                       initfreqs="numeric",
                       models="list",  # named list of 'contextModel's
                       modelnames="character", # list of model names named by tip, node names of tree
                       likfun="function",
                       results="list",
                       invocation="character"
                   )
         )

#' The results of a contextual MCMC on a tree
#'
#' @name contextTreeMCMC-class
#' @rdname contextTreeMCMC-class
#' @exportClass contextTreeMCMC
setClass("contextTreeMCMC", 
        representation(
                       mutrates.prior.mean="numeric",
                       mutrates.prior.sd="numeric",
                       selcoef.prior.mean="numeric",
                       selcoef.prior.sd="numeric",
                       fixfn.params.prior.mean="numeric"
                       fixfn.params.prior.sd="numeric"
                   ),
             contains="contextTree"
         )

#' @rdname context-methods
#' @aliases dimnames,context-method
#' @export dimnames
setMethod("dimnames", signature=c(x="context"), definition=function (x) { dimnames(x@counts) } )
#' Methods for contextTree
#'
#' @rdname contextTree-methods
#' @aliases dimnames,contextTree-method
#' @export dimnames
setMethod("dimnames", signature=c(x="contextTree"), definition=function (x) { dimnames(x@counts) } )
#' @rdname context-methods
#' @aliases counts,context-method
#' @export counts
setMethod("counts", signature=c(x="context"), definition=function (x) { counts(x@counts) } )
#' @rdname contextTree-methods
#' @aliases counts,contextTree-method
#' @export counts
setMethod("counts", signature=c(x="contextTree"), definition=function (x) { counts(x@counts) } )
#' @rdname context-methods
#' @aliases countframe,context-method
#' @export countframe
setMethod("countframe", signature=c(x="context"), definition=function (x) countframe(x@counts))
#' @rdname contextTree-methods
#' @aliases countframe,contextTree-method
#' @export countframe
setMethod("countframe", signature=c(x="contextTree"), definition=function (x) countframe(x@counts))


# extract window lengths from these objects:
#' Methods to get window lengths.
#'
#' @name longwin
#' @rdname context-methods
#' @exportMethod longwin
setGeneric("longwin", function(x) { standardGeneric("longwin") })
#' @name shortwin
#' @rdname context-methods
#' @exportMethod shortwin
setGeneric("shortwin", function(x) { standardGeneric("shortwin") })
#' @name leftwin
#' @rdname context-methods
#' @exportMethod leftwin
setGeneric("leftwin", function(x) { standardGeneric("leftwin") })
#' @name rightwin
#' @rdname context-methods
#' @exportMethod rightwin
setGeneric("rightwin", function(x) { longwin(x)-shortwin(x)-leftwin(x) })

#' Methods for genmatix objects
#'
#' @rdname genmatrix-methods
#' @aliases longwin,genmatrix-method
#' @export longwin
setMethod("longwin", signature=c(x="genmatrix"), definition=function(x) { nchar(rownames(x)[1]) } )
#' @rdname tuplecounts-methods
#' @aliases longwin,tuplecounts-method
#' @export longwin
setMethod("longwin", signature=c(x="tuplecounts"), definition=function(x) { nchar(rownames(x@counts)[1]) } )
#' @rdname context-methods
#' @aliases longwin,context-method
#' @export longwin
setMethod("longwin", signature=c(x="context"), definition=function(x) { longwin(x@counts) } )
#' @rdname contextTree-methods
#' @aliases longwin,contextTree-method
#' @export longwin
setMethod("longwin", signature=c(x="contextTree"), definition=function(x) { longwin(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases shortwin,tuplecounts-method
#' @export shortwin
setMethod("shortwin", signature=c(x="tuplecounts"), definition=function(x) { unique(sapply(lapply(lapply(colpatterns(x),"[",1),as.character),nchar)) } )
#' @rdname context-methods
#' @aliases shortwin,context-method
#' @export shortwin
setMethod("shortwin", signature=c(x="context"), definition=function(x) { shortwin(x@counts) } )
#' @rdname contextTree-methods
#' @aliases shortwin,contextTree-method
#' @export shortwin
setMethod("shortwin", signature=c(x="contextTree"), definition=function(x) { shortwin(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases leftwin,tuplecounts-method
#' @export leftwin
setMethod("leftwin", signature=c(x="tuplecounts"), definition=function(x) { x@leftwin } )
#' @rdname context-methods
#' @aliases leftwin,context-method
#' @export leftwin
setMethod("leftwin", signature=c(x="context"), definition=function(x) { leftwin(x@counts) } )
#' @rdname contextTree-methods
#' @aliases leftwin,contextTree-method
#' @export leftwin
setMethod("leftwin", signature=c(x="contextTree"), definition=function(x) { leftwin(x@counts) } )

##############
# convenience functions

#' @name nmuts
#' @rdname context-methods
#' @exportMethod nmuts
setGeneric("nmuts", function (x) { standardGeneric("nmuts") } )
#' @name nsel
#' @rdname context-methods
#' @exportMethod nsel
setGeneric("nsel", function (x) { standardGeneric("nsel") } )
#' @name fixparams
#' @rdname context-methods
#' @exportMethod fixparams
setGeneric("fixparams", function (x) { standardGeneric("fixparams") } )
#' @rdname genmatrix-methods
#' @aliases nmuts,genmatrix-method
#' @export nmuts
setMethod("nmuts", signature=c(x="genmatrix"), definition=function (x) { length(x@mutpats) } )
#' @rdname genmatrix-methods
#' @aliases nsel,genmatrix-method
#' @export nsel
setMethod("nsel", signature=c(x="genmatrix"), definition=function (x) { length(x@selpats) } )
#' @rdname genmatrix-methods
#' @aliases fixparams,genmatrix-method
#' @export fixparams
setMethod("fixparams", signature=c(x="genmatrix"), definition=function (x) { (setdiff(names(as.list(formals(x@fixfn))),"..."))[-1] } )
#' Methods for contextModel objects
#'
#' @rdname contextModel-methods
#' @aliases nmuts,contextModel-method
#' @export nmuts
setMethod("nmuts", signature=c(x="contextModel"), definition=function (x) { nmuts(x@genmatrix) } )
#' @rdname contextModel-methods
#' @aliases nsel,contextModel-method
#' @export nsel
setMethod("nsel", signature=c(x="contextModel"), definition=function (x) { nsel(x@genmatrix) } )
#' @rdname contextModel-methods
#' @aliases fixparams,contextModel-method
#' @export fixparams
setMethod("fixparams", signature=c(x="contextModel"), definition=function (x) { fixparams(x@genmatrix) } )
#' @rdname contextTree-methods
#' @aliases nmuts,contextTree-method
#' @export nmuts
setMethod("nmuts", signature=c(x="contextTree"), definition=function (x) { sapply( x@models, nmuts ) } )
#' @rdname contextTree-methods
#' @aliases nsel,contextTree-method
#' @export nsel
setMethod("nsel", signature=c(x="contextTree"), definition=function (x) { sapply( x@models, nsel ) } )
#' @rdname contextTree-methods
#' @aliases fixparams,contextTree-method
#' @export fixparams
setMethod("fixparams", signature=c(x="contextTree"), definition=function (x) { lapply( x@models, fixparams ) } )

#######
# methods related to model fitting

#' @rdname context-methods
#' @aliases coef,context-method
#' @export coef
setMethod("coef", signature=c(object="context"), definition=function (object) {
          coef <- c( object@mutrates, object@selcoef, object@params )
          names(coef) <- c( mutnames( object@genmatrix@mutpats, pretty=TRUE ), selnames( object@genmatrix@selpats, pretty=TRUE ), names(object@params) )
          return(coef) } )

#' @rdname contextModel-methods
#' @aliases coef,contextModel-method
#' @export coef
setMethod("coef", signature=c(object="contextModel"), definition=function (object) {
          coef <- c( object@mutrates, object@selcoef, object@params )
          names(coef) <- c( mutnames( object@genmatrix@mutpats, pretty=TRUE ), selnames( object@genmatrix@selpats, pretty=TRUE ), names(object@params) )
          return(coef) } )

#' @rdname contextTree-methods
#' @aliases coef,contextTree-method
#' @export coef
setMethod("coef", signature=c(object="contextTree"), definition=function (object) { 
            tlens <- object@tree$edge.length
            names(tlens) <- paste("tlen",edge.labels(object@tree),sep='.')
            c( tlens, do.call( c, lapply( object@models, coef ) ) ) 
        } )
#' @rdname tuplecounts-methods
#' @aliases rowSums,tuplecounts-method
#' @export rowSums
setMethod("rowSums", signature=c(x="tuplecounts"), definition=function (x) { Matrix::rowSums(x@counts) } )
#' @rdname tuplecounts-methods
#' @aliases image,tuplecounts-method
#' @export image
setMethod("image", signature=c(x="tuplecounts"), definition=function (x) { image(x@counts) } )
#' @rdname context-methods
#' @aliases rowSums,context-method
#' @export rowSums
setMethod("rowSums", signature=c(x="context"), definition=function (x) { Matrix::rowSums(x@counts@counts) } )
#' @rdname contextTree-methods
#' @aliases rowSums,contextTree-method
#' @export rowSums
setMethod("rowSums", signature=c(x="contextTree"), definition=function (x) { Matrix::rowSums(x@counts@counts) } )
#' @rdname context-methods
#' @aliases fitted,context-method
#' @export fitted
setMethod("fitted", signature=c(object="context"), definition=function (object,...) { predictcounts.context(object,...) } )
#' @rdname contextTree-methods
#' @aliases fitted,contextTree-method
#' @export fitted
setMethod("fitted", signature=c(object="contextTree"), definition=function (object,...) { 
              predictcounts.contextTree(object, ...)
        } )
#' @rdname context-methods
#' @aliases residuals,context-method
#' @export residuals
setMethod("residuals", signature=c(object="context"), definition=function (object,...) { resid.context(object,...) } )


## This makes everything go all to hell, for some reason:
# setMethod("-",signature=c("tuplecounts","tuplecounts"), definition=function(e1,e2) { e1@counts <- (e1@counts-e2@counts); return(e1) } )
# setMethod("-",signature=c("tuplecounts","ANY"), definition=function(e1,e2) { new("tuplecounts",leftwin=e1@leftwin,counts=Matrix(e1@counts-e2),bases=e1@bases) } )
# setMethod("-",signature=c("ANY","tuplecounts"), definition=function(e1,e2) { new("tuplecounts",leftwin=e2@leftwin,counts=Matrix(e1-e2@counts),bases=e2@bases) } )

