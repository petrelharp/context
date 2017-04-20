library(simcontext)
library(contextual)
library(testthat)

context("An extreme model: has X's as VERY deleterious; O's as VERY beneficial and  sp2 on a LONG branch")

varying_selcoef_config <- '
{
    "comment" : "Test assignment of weights to selection coefficients.",
    "tree" : [ "(sp1 : 1.0, sp2 : 20.0) root;" ],
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 0.5, 0.5 ],
    "model" : {
        "mutpats" : [
            [ [ "O", "X" ] ],
            [ [ "X", "O" ] ]
        ],
        "mutrates" : [ 0.5, 0.5 ],
        "comment": "note unequal weights in first line",
        "selpats" : [
            { "O" : 20.0, "X" : -20.0 }
        ],
        "selcoef" : [ 1.0 ],
        "fixfn" : "ising.fixfn",
        "fixfn.params" : []
    },
    "sp1" : "model",
    "sp2" : "model"
}
'

config <- parse.models( treeify.config( read.config(json=varying_selcoef_config) ) )

simseqs <- simseq.tree(100,config)

longpats <- getpatterns(4,bases=config$bases)
shortpats <- getpatterns(1,bases=config$bases)

# counts
shorts.1 <- counttrans(shortpats, shortpats, simseqs=simseqs$sp1 )
shorts.2 <- counttrans(shortpats, shortpats, simseqs=simseqs$sp2 )

test_that("for sp1, should have no O->X",
    expect_equal( counts(shorts.1)["O","X"], 0 )
)
test_that("for sp2, should have no X's at all", {
    expect_equal(counts(shorts.2)["O","X"], 0 )
    expect_equal(counts(shorts.2)["X","X"], 0 )
} )



