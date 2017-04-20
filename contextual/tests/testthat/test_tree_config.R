
library(contextual)
library(testthat)

context("Check expected numbers under an extreme model.")

## Pattern should be:
# root: X
#  root -> an1: X -> X
#  root -> an2: X -> O
#   an1 -> sp1: X -> O
#   an1 -> sp2: X -> O
#   an2 -> sp3: O -> O
#   an2 -> sp4: O -> X

another_big_tree <- '
{
    "comment" : "Complicated tree, for testing.",
    "tree" : [ "( (sp1 : 0.1, sp2 : 0.1) an1 : 1.0, (sp3 : 1.0, sp4 : 0.2) an2 : 3.0 ) root;" ],
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 1.0, 0.0 ],
    "mutpats" : [
        [ [ "X", "O" ] ]
    ],
    "mutrates" : [ 10 ],
    "selpats" : [ [ "X" ] ],
    "selcoef" : [ 0.01 ],
    "fixfn" : "popgen.fixfn",
    "fixfn.params" : { "Ne" : 100 },
    "an1" : { 
            "comment" : "do nothing",
            "selfactors" : [ 100 ]
        },
    "an2" : { 
            "comment" : "default, X -> O"
        },
    "sp1" : { 
            "comment" : "X -> O",
            "selpats" : [ [ "O" ] ],
            "selfactors" : [ [ 5 ] ]
        },
    "sp2" : { 
            "comment" : "X -> O",
            "mutpats" : [
                [ [ "X", "O" ] ],
                [ [ "O", "X" ] ]
            ],
            "mutrates" : [ 1, 1 ],
            "selpats" : [ [ "O" ] ],
            "selcoef" : [ 1.0 ]
        },
    "sp3" : { 
            "comment" : "do nothing",
            "mutpats" : [ [ "O", "X" ] ],
            "mutrates" : [ 10 ],
            "selpats" : [ [ "O" ] ],
            "selcoef" : [ 100 ]
        },
    "sp4" : {
            "comment" : "O -> X",
            "mutpats" : [ [ "O", "X" ] ],
            "mutrates" : [ 30 ]
        }
}
'

config <- parse.models( treeify.config( read.config(json=another_big_tree) ) )

clades <- c("an1","an2","sp1","sp2","sp3","sp4")
names(clades) <- clades

genmat.list <- lapply( clades, function (clade) {
        with( list2env(config[[clade]]), 
                makegenmatrix(mutpats,selpats,patlen=1,patterns=bases,fixfn=fixfn,mutrates=mutrates,selcoef=selcoef,selfactors=selfactors,bases=bases,Ne=fixfn.params$Ne)
            )
    } )


pf <- function (ds) { do.call(config$fixfn,c(list(ds),config$fixfn.params)) }

expected.genmats <- cbind( # X->X O->X X->O O->O
                    an1 = c(  0,  0,   0,   0),
                    an2 = c(  0,  0, 10*pf(-.01), 0),
                    sp1 = c(  0,  0, 10*pf(.05),   0),
                    sp2 = c(  0,  0,    pf(1),   0),
                    sp3 = c(  0,  0,   0,   0),
                    sp4 = c(  0, 30*pf(.01),   0,   0)
        )

test_that("Counts match expected numbers", {
    expect_true( all(which(expected.genmats>0)==which(sapply(genmat.list,as.vector)>1e-8)) )
    expect_equal(
            as.vector(expected.genmats),
            as.vector(sapply(genmat.list,as.vector))
    )
})
