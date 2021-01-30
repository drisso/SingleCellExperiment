# This tests the applySCE function.
# library(testthat); library(SingleCellExperiment); source("setup.R"); source("test-applySCE.R")

TESTFUN <- function(x, a=1, b=2, c=3) {
    list(X=x, ARGS=list(A=a, B=b, C=c))
}

identity2 <- function(x) {
    if (is(x, "SingleCellExperiment")) removeAltExps(x) else x
}

test_that("applySCE works as expected", {
    def <- applySCE(loaded, FUN=TESTFUN)
    expect_identical(def[[1]]$X, loaded)
    expect_identical(def[[2]]$X, altExp(loaded, "Spike"))
    expect_identical(def[[3]]$X, altExp(loaded, "Protein"))
    expect_identical(unique(lapply(def, function(x) x$ARGS)), list(list(A=1, B=2, C=3)))
    expect_identical(names(def), c("", "Spike", "Protein"))

    # Trying with custom args.
    out <- applySCE(loaded, FUN=TESTFUN, MAIN.ARGS=list(a=-1),
        ALT.ARGS=list(Spike=list(b=-2, c=-10)), a=-100)

    expect_identical(out[[1]]$ARGS, list(A=-1, B=2, C=3))
    expect_identical(out[[2]]$ARGS, list(A=-100, B=-2, C=-10))
    expect_identical(out[[3]]$ARGS, list(A=-100, B=2, C=3))

    # Other choices of WHICH work.
    def <- applySCE(loaded, FUN=TESTFUN, WHICH="Spike")
    short <- applySCE(loaded, FUN=TESTFUN, WHICH=1)
    expect_identical(def, short)

    # Check that the right names are pulled by simplifyToSCE.
    def <- applySCE(loaded, FUN=identity2)
    short <- applySCE(loaded, FUN=identity2, WHICH=1:2)
    expect_identical(def, short)
})

test_that("applySCE errors out nicely", {
    expect_error(out <- applySCE(loaded, FUN=function(x) stop("YAY")), "YAY")
    expect_error(out <- applySCE(loaded, FUN=function(x) stop("YAY")), "failed on the main Experiment", fixed=TRUE)
    expect_error(out <- applySCE(loaded, MAIN.ARGS=NULL, FUN=function(x) stop("YAY")), "failed on alternative Experiment", fixed=TRUE)
})

test_that("simplification works correctly", {
    results <- applySCE(loaded, FUN=identity2)
    expect_identical(results, loaded)

    # Manual simplification works.
    raw <- applySCE(loaded, FUN=identity2, WHICH=c("Spike", "Protein"), SIMPLIFY=FALSE)
    expect_identical(raw, list(removeAltExps(loaded), Spike=altExp(loaded), Protein=altExp(loaded, 2)))
    expect_identical(results, simplifyToSCE(raw)) # inference works correctly.

    # Simplifies even when main is absent.
    raw <- applySCE(loaded, FUN=identity2, MAIN.ARGS=NULL)
    expect_identical(nrow(raw), 0L)
    expect_identical(altExps(raw), altExps(loaded)) 
})

test_that("simplification fails correctly", {
    raw <- applySCE(loaded, FUN=identity2, SIMPLIFY=FALSE, WHICH=c("Spike", "Spike"))
    expect_warning(results <- simplifyToSCE(raw), "multiple references.*Spike")
    expect_error(results <- simplifyToSCE(raw, warn.level=3), "multiple references.*Spike")
    expect_null(results)

    raw <- applySCE(loaded, FUN=identity2, SIMPLIFY=FALSE)
    expect_error(simplifyToSCE(raw, which.main=1:2), "length(which.main)", fixed=TRUE)
    expect_error(simplifyToSCE(unname(raw)), "multiple")

    raw <- applySCE(loaded, FUN=function(x) 1, SIMPLIFY=FALSE, WHICH=1)
    expect_warning(results <- simplifyToSCE(raw, warn.level=1), NA)
    expect_warning(results <- simplifyToSCE(raw), "not a SummarizedExperiment")
    expect_null(results)

    set.seed(1000)
    raw <- applySCE(loaded, FUN=function(x) { x[,sample(ncol(x), sample(ncol(x), 1))] }, SIMPLIFY=FALSE)
    expect_warning(results <- simplifyToSCE(raw, warn.level=0), NA)
    expect_warning(results <- simplifyToSCE(raw), "columns are different")
    expect_null(results)

    raw <- applySCE(loaded, FUN=function(x) { as(x, "SummarizedExperiment") }, SIMPLIFY=FALSE, WHICH=1)
    expect_warning(results <- simplifyToSCE(raw, warn.level=1), NA)
    expect_warning(results <- simplifyToSCE(raw), "not a SingleCellExperiment")
    expect_null(results)
})
