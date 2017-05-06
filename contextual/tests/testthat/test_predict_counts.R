library(contextual)
library(testthat)

config <- parse.models(list(
    "tree" = c( "(sp1 : 1.8, (sp2 : 1.2, sp3 : 1.0) an1 : 0.5 ) root;" ), 
    "bases" = c( "X", "O" ),
    "initfreqs" = c( 0.5, 0.5),
    "sp1" = list(
        "mutpats" = list(
            list( c("X", "O" ) ),
            list( c("O", "X" ) ),
            list( c("XX", "XO" ), c( "OO", "XO" ) )
        ),
        "mutrates" = c( 10.0, 0.01, 0.01 )
    ),
    "sp2" = "sp1",
    "sp3" = "sp1",
    "an1" = "sp1"
))

# nontree

counts <- new("tuplecounts", leftwin = 0, 
              counts = new("dgCMatrix", 
                           i = c(0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L), 
                           p = c(0L, 4L, 8L, 12L, 16L), Dim = c(4L, 4L), 
                           Dimnames = list(c("XX", "OX", "XO", "OO"), c("XX", "OX", "XO", "OO")), 
                           x = c(210, 59, 50, 16, 24, 141, 13, 41, 21, 11, 153, 35, 4, 35, 30, 156), 
                           factors = list()), 
              bases = c("O", "X"), rowtaxon = "long", 
              colpatterns = structure(list(short = structure(c(4L, 2L, 3L, 1L), 
                                                             .Label = c("OO", "OX", "XO", "XX"), 
                                                             class = "factor")), 
                                      .Names = "short", row.names = c(NA, -4L), class = "data.frame"))


G <- do.call(makegenmatrix, c(list(patlen=2), config$sp1))
projmatrix <- collapsepatmatrix(ipatterns=getpatterns(2,config$bases),
                                leftwin=0, shortwin=1, rightwin=1, bases=config$bases)

model <- new("context",
               genmatrix=G,
               projmatrix=projmatrix,
               mutrates=config$sp1$mutrates,
               selcoef=numeric(0),
               params=numeric(0),
               counts=counts,
               likfun=function(x) {1},
               results=list(),
               invocation=""
            )

pred <- fitted(model)

test_that("Predicted counts are sensible.", {
              expect_equal(sum(pred@counts), sum(model@counts@counts))
              expect_true(all(pred@counts[,"X"]<1.0))
        } )


# tree

counts <- new("tuplecounts", 
              leftwin = 0, 
              counts = new("dgCMatrix", 
                       i = c(0L, 1L, 2L, 0L, 3L, 0L, 1L, 2L, 3L, 1L, 0L, 1L, 0L, 1L, 2L, 3L, 0L, 1L, 3L, 0L, 1L, 2L, 3L, 0L, 3L, 0L, 1L, 2L, 3L, 0L, 2L, 3L, 2L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 1L, 2L, 3L) , p = c(0L, 3L, 5L, 9L, 10L, 12L, 16L, 16L, 19L, 23L, 25L, 29L, 32L, 33L, 37L, 41L, 44L), 
                       Dim = c(4L, 16L), 
                       Dimnames = list(c("XX", "OX", "XO", "OO"), c("XX.XX", "OX.XX", "XO.XX", "OO.XX", "XX.OX", "OX.OX", "XO.OX", "OO.OX", "XX.XO", "OX.XO", "XO.XO", "OO.XO", "XX.OO", "OX.OO", "XO.OO", "OO.OO")), 
                       x = c(4, 4, 5, 1, 1, 2, 1, 1, 1, 1, 4, 1, 1, 8, 3, 6, 1, 1, 2, 1, 2, 1, 1, 2, 1, 3, 3, 7, 3, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 3, 1, 4, 4), 
                       factors = list()), 
              bases = c("X", "O"), 
              rowtaxon = "sp1", 
              colpatterns = structure(list(
                       sp2 = c("XX", "OX", "XO", "OO", "XX", "OX", "XO", "OO", "XX", "OX", "XO", "OO", "XX", "OX", "XO", "OO"), 
                       sp3 = c("XX", "XX", "XX", "XX", "OX", "OX", "OX", "OX", "XO", "XO", "XO", "XO", "OO", "OO", "OO", "OO")), 
                      .Names = c("sp2", "sp3"), 
                      out.attrs = structure(list( dim = c(4L, 4L), dimnames = structure(list(Var1 = c("Var1=XX", "Var1=OX", "Var1=XO", "Var1=OO"), Var2 = c("Var2=XX", "Var2=OX", "Var2=XO", "Var2=OO")), .Names = c("Var1", "Var2"))), .Names = c("dim", "dimnames")), 
                  class = "data.frame", row.names = c(NA, -16L))
          )

projmatrix.2 <- collapsepatmatrix(ipatterns=getpatterns(2,config$bases),
                                leftwin=0, shortwin=2, rightwin=0, bases=config$bases)

model.2 <- new("context",
               genmatrix=G,
               projmatrix=projmatrix.2,
               mutrates=config$sp1$mutrates,
               selcoef=numeric(0),
               params=numeric(0),
               counts=counts,
               likfun=function(x) {1},
               results=list(),
               invocation=""
            )

treemodel <- new("contextTree",
                 counts=counts,
                 tree=config$tree,
                 initfreqs=c(1/2,1/2),
                 models=list(sp1=model.2),
                 modelnames=config.dereference(config, 
                                               nodenames(config$tree)),
                 likfun=function(x) {1},
                 results=list(),
                 invocation=""
            )

# conditioned on sp1 counts
pred.sp1 <- fitted(treemodel, rowtaxon="sp1", coltaxa=c("sp2","sp3"))
# not conditioned on these
pred.sp1.mean <- fitted(treemodel, rowtaxon="sp1", coltaxa=c("sp2","sp3"), 
                        initcounts=sum(counts@counts))

test_that("Predicted tree counts are sensible.", {
              expect_equal(sum(pred.sp1@counts), sum(treemodel@counts@counts))
              expect_equal(sum(pred.sp1.mean@counts), sum(treemodel@counts@counts))
              expect_true(all(pred.sp1@counts[,-match("OO,OO",colnames(pred.sp1@counts))]<1.0))
              expect_true(pred.sp1.mean@counts["OO","OO,OO"] > sum(treemodel@counts@counts)-1)
        } )
