library(contextual)
library(testthat)

context("Reading in configuration files")

model_json <- '
{
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 0.5, 0.5 ],
    "mutpats" : [
        [ [ "XO", "OX" ] ]
    ],
    "mutrates" : [ 1 ],
    "mutrates.scale" : [ 0.1 ]
}
'
model_config <- read.config(json=model_json)


expect_list_equal <- function (x,y) {
    expect_equal(length(x), length(y))
    expect_equal(names(x), names(y))
    for (n in names(x)) {
        expect_equal(x[[n]], y[[n]])
    }
}

expect_list_equal(model_config,
                    list(bases = c("X", "O"), 
                         initfreqs = c(0.5, 0.5), 
                         mutpats = list(list(c("XO", "OX"))), 
                         mutrates = 1L, 
                         mutrates.scale = 0.1))

config <- treeify.config(model_config,tlen=10)

expect_list_equal(config,
                    list(tree = structure(list(edge = structure(c(2L, 1L), .Dim = 1:2), Nnode = 1L, tip.label = "tip", node.label = "root", edge.length = 10), 
                                          .Names = c("edge", "Nnode", "tip.label", "node.label", "edge.length"), class = "phylo"), 
                         bases = c("X", "O"), 
                         tip = structure(list(bases = c("X", "O"), initfreqs = c(0.5, 0.5), mutpats = list(list(c("XO", "OX"))), mutrates = 1L, mutrates.scale = 0.1), 
                                         .Names = c("bases", "initfreqs", "mutpats", "mutrates", "mutrates.scale")), 
                         initfreqs = c(0.5, 0.5)) )

config <- parse.models(config)

expect_list_equal(config, 
                    list(tree = structure(list(edge = structure(c(2L, 1L), .Dim = 1:2), Nnode = 1L, tip.label = "tip", node.label = "root", edge.length = 10), 
                                          .Names = c("edge", "Nnode", "tip.label", "node.label", "edge.length"), class = "phylo"), 
                         bases = c("X", "O"), 
                         tip = structure(list(bases = c("X", "O"), initfreqs = c(0.5, 0.5), mutpats = list(list(c("XO", "OX"))), 
                                              mutrates = 1L, mutrates.scale = 0.1, selpats = list(), fixfn.params = list(), fixfn.params.scale = list(), 
                                              selfactors = list(), fixfn = function (...) { 1 }, selcoef = numeric(0), selcoef.scale = numeric(0)), 
                                         .Names = c("bases", "initfreqs", "mutpats", "mutrates", "mutrates.scale", "selpats", "fixfn.params", 
                                                    "fixfn.params.scale", "selfactors", "fixfn", "selcoef", "selcoef.scale")), 
                         initfreqs = c(0.5, 0.5), 
                         .models = structure("tip", .Names = "tip")) )


context("Reading in tree configuration")

tree_model_json <- '
{
    "tree" : [ "(sp1 : 0.8, sp2 : 1.2) root;" ],
    "bases" : [ "X", "O" ],
    "initfreqs" : [ 0.5, 0.5 ],
    "sp1" : {
        "comment" : "Each X moves right.",
        "mutpats" : [
            [ [ "XO", "OX" ] ]
        ],
        "mutrates" : [ 1 ]
    },
    "sp2" : {
        "comment" : "Each X moves left.",
        "mutpats" : [
            [ [ "OX", "XO" ] ]
        ],
        "mutrates" : [ 1 ]
    }
}
'
tree_model_config <- read.config(json=tree_model_json)

expect_list_equal( tree_model_config,
                    list(tree = "(sp1 : 0.8, sp2 : 1.2) root;", 
                         bases = c("X", "O"), 
                         initfreqs = c(0.5, 0.5), 
                         sp1 = structure(list(comment = "Each X moves right.", mutpats = list(list(c("XO", "OX"))), mutrates = 1L), 
                                         .Names = c("comment", "mutpats", "mutrates")), 
                         sp2 = structure(list(comment = "Each X moves left.", mutpats = list(list(c("OX", "XO"))), mutrates = 1L), 
                                         .Names = c("comment", "mutpats", "mutrates"))) )

tree_config <- treeify.config(tree_model_config)

expect_list_equal( tree_config, 
                    list(tree = structure(list(edge = structure(c(3L, 3L, 1L, 2L), .Dim = c(2L, 2L)), Nnode = 1L, tip.label = c("sp1", "sp2"), 
                                               edge.length = c(0.8, 1.2), node.label = "root"), 
                                          .Names = c("edge", "Nnode", "tip.label", "edge.length", "node.label"), class = "phylo", order = "cladewise"), 
                        bases = c("X", "O"), 
                        initfreqs = c(0.5, 0.5), 
                        sp1 = structure(list( comment = "Each X moves right.", mutpats = list(list( c("XO", "OX"))), mutrates = 1L), 
                                        .Names = c("comment", "mutpats", "mutrates")), 
                        sp2 = structure(list(comment = "Each X moves left.", mutpats = list(list(c("OX", "XO"))), mutrates = 1L), 
                                        .Names = c("comment", "mutpats", "mutrates"))) )

tree_config <- parse.models(tree_config)

expect_list_equal( tree_config, 
                    list(tree = structure(list(edge = structure(c(3L, 3L, 1L, 2L), .Dim = c(2L, 2L)), Nnode = 1L, tip.label = c("sp1", "sp2"), 
                                               edge.length = c(0.8, 1.2), node.label = "root"), 
                                          .Names = c("edge", "Nnode", "tip.label", "edge.length", "node.label"), class = "phylo", order = "cladewise"), 
                        bases = c("X", "O"), 
                        initfreqs = c(0.5, 0.5), 
                        sp1 = structure(list( comment = "Each X moves right.", mutpats = list(list( c("XO", "OX"))), mutrates = 1L, 
                                             selpats = list(), bases = c("X", "O"), fixfn.params = list(), fixfn.params.scale = list(), 
                                             selfactors = list(), fixfn = function (...) { 1 }, selcoef = numeric(0), selcoef.scale = numeric(0)), 
                                        .Names = c("comment", "mutpats", "mutrates", "selpats", "bases", "fixfn.params", "fixfn.params.scale", 
                                                   "selfactors", "fixfn", "selcoef", "selcoef.scale")), 
                        sp2 = structure(list(comment = "Each X moves left.", mutpats = list(list(c("OX", "XO"))), mutrates = 1L, 
                                             selpats = list(), bases = c("X", "O"), fixfn.params = list(), fixfn.params.scale = list(), 
                                             selfactors = list(), fixfn = function (...) { 1 }, selcoef = numeric(0), selcoef.scale = numeric(0)), 
                                        .Names = c("comment", "mutpats", "mutrates", "selpats", "bases", "fixfn.params", "fixfn.params.scale", 
                                                   "selfactors", "fixfn", "selcoef", "selcoef.scale")), 
                        .models = structure(c("sp1", "sp2"), .Names = c("sp1", "sp2"))) )

####
context("selection patterns")

model_json <- '
{
    "bases" : ["X","O"],
    "mutpats" : [ [ ["XO","OX"] ], [ ["OX", "XO"] ] ],
    "mutrates" : [ 2, 1 ],
    "selpats" : { "double" : { "XX" : [2], "OO" : [1] }, "single" : { "X" : [1] } },
    "selcoef" : { "double" : 5, "single" : 1 }
}
'

model_config <- read.config(json=model_json)

expect_list_equal( model_config,
                    list(bases = c("X", "O"), 
                         mutpats = list(list(c("XO", "OX")), list(c("OX", "XO"))), 
                         mutrates = c(2L, 1L), 
                         selpats = list( double = c("XX", "OO"), single = "X"), 
                         selcoef = list(double = 5L, single = 1L), 
                         selfactors = list(double = c(XX=2, OO=1), single = c(X=1)) ))

config <- parse.models(treeify.config(model_config,tlen=10))

test_that("selfactors are correctly parsed", {
        expect_list_equal( config,
                    list(
                         tree = structure(list(edge = structure(c(2L, 1L), .Dim = 1:2), Nnode = 1L, tip.label = "tip", node.label = "root", 
                                               edge.length = 10), .Names = c("edge", "Nnode", "tip.label", "node.label", "edge.length"), class = "phylo"), 
                         bases = c("X", "O"), 
                         tip = list(
                                  bases = c("X", "O"), 
                                  mutpats = list( list(c("XO", "OX")), list(c("OX", "XO"))), 
                                  mutrates = c(2L, 1L), 
                                  selpats = structure(list(double = c("XX", "OO"), single = "X"), .Names = c("double", "single")), 
                                  selcoef = structure(list(double = 5L, single = 1L), .Names = c("double", "single")), 
                                  selfactors = structure(list(double = structure(c(2, 1), .Names = c("XX", "OO")), single = structure(1, .Names = "X")), .Names = c("double", "single")), 
                                  fixfn.params = list(), fixfn.params.scale = list(), 
                                  fixfn = function (...) { 1 }),
                  initfreqs = NULL, 
                  .models = structure("tip", .Names = "tip")) )
})
