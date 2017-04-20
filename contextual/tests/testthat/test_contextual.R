library(testthat)

library(contextual)

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


test_that("getmutpats returns stored", {
    expect_equal(
        getmutpats.out,
        getmutpats(bases=bases,patlen=2)
    )
})

test_that("mutnames works", {
    expect_equal(mutnames(mutpats),"XO->OX")
})


# if functions have been updated
# save(mutpatchanges.out, getmutpats.out, file=tests_data_fname)
