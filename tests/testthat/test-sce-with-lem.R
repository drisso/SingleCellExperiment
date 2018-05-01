# Checks that LEM can be used inside a SCE object.
# library(SingleCellExperiment); library(testthat); source("test-sce-with-lem.R")
context("Test LEM in reducedDim()")

set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
lem <- LinearEmbeddingMatrix(factors, loadings)

v <- matrix(rnorm(20000), ncol=ncells)
u <- matrix(rpois(20000, 5), ncol=ncells)

test_that("SCE construction works correctly with LEM", {

    sce1 <- SingleCellExperiment(assay=SimpleList(counts=u, exprs=v),
                                 reducedDims = SimpleList(rd1 = lem))
    sce2 <- SingleCellExperiment(assay=SimpleList(counts=u, exprs=v),
                                 reducedDims = SimpleList(rd1 = lem, rd2 = factors))


    expect_identical(reducedDim(sce1, "rd1"), reducedDim(sce2, "rd1"))
    expect_identical(reducedDim(sce1, "rd1"), lem)

})

test_that("reduced dimension getters/setters are functioning", {

    sce <- SingleCellExperiment(assay=SimpleList(counts=u, exprs=v))
    reducedDim(sce, "PCA") <- lem
    expect_identical(reducedDim(sce, "PCA"), lem)
    expect_identical(reducedDims(sce), SimpleList(PCA=lem))
    expect_identical(reducedDimNames(sce), "PCA")

})

