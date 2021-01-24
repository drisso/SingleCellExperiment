# This tests the sceApply function.
# library(testthat); library(SingleCellExperiment); source("setup.R"); source("test-sceApply.R")

TESTFUN <- function(x, a=1, b=2, c=3) {
    list(X=x, ARGS=list(A=a, B=b, C=c))
}

test_that("sceApply works as expected", {
    out <- sceApply(loaded, FUN=TESTFUN, which=list(
        AssayInput(assay="logcounts", a=-1),
        ReducedDimInput(type=2, b=-2, c=-10),
        AltExpInput(experiment="Protein")
    ), a=-100)

    expect_identical(out[[1]]$X, assay(loaded, "logcounts"))
    expect_identical(out[[2]]$X, reducedDim(loaded, 2))
    expect_identical(out[[3]]$X, altExp(loaded, "Protein"))

    # Overriding done correctly.
    expect_identical(out[[1]]$ARGS, list(A=-1, B=2, C=3))
    expect_identical(out[[2]]$ARGS, list(A=-100, B=-2, C=-10))
    expect_identical(out[[3]]$ARGS, list(A=-100, B=2, C=3))
})

test_that("sceApply retains names", {
    out <- sceApply(loaded, FUN=TESTFUN, which=list(
        helen=AssayInput(assay="logcounts", a=-1),
        sarah=ReducedDimInput(type=2, b=-2, c=-10),
        vanessa=AltExpInput(experiment="Protein")
    ), a=-100)

    expect_identical(names(out), c("helen", "sarah", "vanessa"))
})
