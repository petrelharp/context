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

# to generate data:
#     library(simcontext)
#     sim <- do.call(simseq, c(list(seqlen=1000, tlen=0.1), config[['sp1']]))
#     simcounts <- counttrans(ipatterns=getpatterns(2,bases=config$bases), 
#                             fpatterns=getpatterns(2,bases=config$bases), simseqs=sim)
#     dput(simcounts)
#     simcounts_long <- counttrans(ipatterns=getpatterns(3,bases=config$bases), 
#                             fpatterns=getpatterns(3,bases=config$bases), simseqs=sim)
#     dput(simcounts_long)

counts <- new("tuplecounts", leftwin = 0, 
              counts = new("dgCMatrix"
                            , i = c(0L, 0L, 1L, 0L, 2L, 0L, 1L, 2L, 3L)
                            , p = c(0L, 1L, 3L, 5L, 9L)
                            , Dim = c(4L, 4L)
                            , Dimnames = list(c("XX", "OX", "XO", "OO"), c("XX", "OX", "XO", "OO"))
                            , x = c(41, 61, 99, 65, 95, 96, 153, 156, 233)
                            , factors = list()),
              bases = c("O", "X"), rowtaxon = "long", 
              colpatterns = structure(list(short = structure(c(4L, 2L, 3L, 1L), 
                                                             .Label = c("OO", "OX", "XO", "XX"), 
                                                             class = "factor")), 
                                      .Names = "short", row.names = c(NA, -4L), class = "data.frame"))

counts_long <- new("tuplecounts"
        , leftwin = 0
        , counts = new("dgCMatrix"
        , i = c(0L, 0L, 1L, 0L, 2L, 0L, 1L, 2L, 3L, 0L, 4L, 0L, 1L, 4L, 5L, 
    0L, 2L, 4L, 6L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L)
        , p = c(0L, 1L, 3L, 5L, 9L, 11L, 15L, 19L, 27L)
        , Dim = c(8L, 8L)
        , Dimnames = list(c("XXX", "OXX", "XOX", "OOX", "XXO", "OXO", "XOO", "OOO"
    ), c("XXX", "OXX", "XOX", "OOX", "XXO", "OXO", "XOO", "OOO"))
        , x = c(12, 8, 21, 15, 24, 17, 29, 29, 46, 11, 18, 26, 28, 27, 50, 
    20, 27, 30, 44, 28, 48, 57, 69, 50, 76, 70, 118)
        , factors = list()
    )
        , bases = c("O", "X")
        , rowtaxon = "long"
        , colpatterns = structure(list(short = c("XXX", "OXX", "XOX", "OOX", "XXO", "OXO", 
    "XOO", "OOO")), .Names = "short", row.names = c(NA, -8L), class = "data.frame")
    )

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

pred <- fitted(model, tlen=1)

test_that("Predicted counts are sensible.", {
              expect_equal(sum(pred@counts), sum(model@counts@counts))
              expect_true(all(pred@counts[,"X"]<1.0))
        } )

##
context("Now predict 3-1-1 counts from a 2-1-0 model.")

longer_G <- do.call(makegenmatrix, c(list(patlen=3), config$sp1))
long_pred <- fitted(model, longwin=3, shortwin=1, leftwin=1, 
                    initcounts=rowSums(counts_long),
                    genmatrix=longer_G)

test_that("Predicting longer counts also sensible.", {
              expect_equal(dim(long_pred), c(2^3, 2))
              expect_equal(sum(long_pred@counts), sum(model@counts@counts)-1)
              expect_equivalent(projectcounts(long_pred, new.leftwin=0, new.longwin=2)@counts, pred@counts)
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
              expect_true(all(pred.sp1@counts[,-match("OO.OO",colnames(pred.sp1@counts))]<1.0))
              expect_true(pred.sp1.mean@counts["OO","OO.OO"] > sum(treemodel@counts@counts)-1)
        } )
