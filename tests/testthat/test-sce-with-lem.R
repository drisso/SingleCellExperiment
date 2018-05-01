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

sce <- SingleCellExperiment(assay=v, reducedDims=SimpleList(PCA=lem)) # A fully loaded object.
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


test_that("subsetting by row works correctly", {

    colData(sce)$indicator <- seq_len(ncells)

    for (j in 1:3) {
        if (j==1L) {
            by.col <- sample(ncells, 20)
            sub.sce <- sce[,by.col]
            expect_identical(colData(sub.sce)$indicator, by.col)
        } else if (j==2L) {
            by.col <- rbinom(ncells, 1, 0.2)==1
            sub.sce <- sce[,by.col]
            expect_identical(colData(sub.sce)$indicator, which(by.col))
        } else if (j==3L) {
            by.col <- colnames(sce)[sample(ncells, 50)]
            sub.sce <- sce[,by.col]
            expect_identical(colData(sub.sce)$indicator, match(by.col, colnames(sce)))
        }
        ind <- colData(sub.sce)$indicator

        expect_identical(reducedDim(sub.sce, "PCA"), lem[ind,,drop=FALSE])

    }

})
