library(testthat)

library("Biobase")
library("Biostrings")
library("IRanges")

source("../context.R")
source("../context-inference-fns.R")

tests_data_fname <- "unit-tests.Rdata"
load(tests_data_fname)

set.seed(42)
bases <- c("X","O")
mutpats <- list(
    list( c("XO","OX") )
    )
mutrates <- 1
selpats <- list()
selcoef <- numeric(0)
fixfn <- function (...) { 1 }


context("Testing basic functions")

test_that("getpatterns works", {
    expect_equal(getpatterns(2,bases), c("XX","OX","XO","OO"))
})

test_that("mutpatchanges returns stored", {
    expect_equal(
        mutpatchanges.out,
        mutpatchanges(
            list(
                list( c("XO","OX"), c("XO","XB") ),
                list( c("A","X") )
                )
        )
    )
})

getmutpats(2)

test_that("getmutpats returns stored", {
    expect_equal(
        getmutpats.out,
        getmutpats(2)
    )
})

test_that("npatterns works", {
    expect_equal(npatterns(2,bases), 4)
})

test_that("mutanmes works", {
    expect_equal(mutnames(mutpats),"XO->OX")
})


# if functions have been updated
# save(mutpatchanges.out, getmutpats.out, file=tests_data_fname)
