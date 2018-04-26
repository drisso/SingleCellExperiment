# Checks for proper construction and get/setting of the slots of a LinearEmbeddingMatrix. 
# library(SingleCellExperiment); library(testthat); source("test-lem-class.R")

set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
lem <- LinearEmbeddingMatrix(factors, loadings)

test_that("LEM construction works correctly", {
    expect_identical(sampleFactors(lem), factors)
    expect_identical(featureLoadings(lem), loadings)
    expect_identical(nrow(factorData(lem)), ncol(lem))
    expect_identical(ncol(factorData(lem)), 0L)

    # Trying again with an initialized factorData.
    fd <- DataFrame(YAY=seq_len(ncol(factors)))
    lem2 <- LinearEmbeddingMatrix(factors, loadings, factorData=fd)
    expect_identical(factorData(lem2), fd)
    expect_identical(lem2$YAY, fd$YAY)

    # Should throw errors if it doesn't make sense.
    expect_error(LinearEmbeddingMatrix(factors, loadings[,1:2]), "must have the same number of columns")
    expect_error(LinearEmbeddingMatrix(factors, loadings, DataFrame(YAY=2)), "one row per factor")
})

test_that("LEM setters work correctly", {
    lem2 <- lem
    sampleFactors(lem2) <- factors * 2
    expect_identical(sampleFactors(lem2), factors * 2)

    lem2 <- lem
    featureLoadings(lem2) <- -loadings
    expect_identical(featureLoadings(lem2), -loadings)

    to.add <- LETTERS[seq_len(ncol(factors))]
    lem2 <- lem
    factorData(lem2)$whee <- to.add
    expect_identical(factorData(lem2)$whee, to.add)

    # Factor Data settters.
    lem2$whee <- NULL
    expect_identical(lem2$whee, NULL)
    lem2$whee <- 42
    expect_identical(lem2$whee, rep(42, ncol(lem2)))

    # Dimnames setters.
    lem2 <- lem
    colnames(lem2) <- to.add
    expect_identical(colnames(lem2), to.add)
    cell.names <- paste0("Cell", seq_len(nrow(factors)))
    rownames(lem2) <- cell.names
    expect_identical(rownames(lem2), cell.names)

    # Should throw errors if it doesn't make sense.
    expect_error(sampleFactors(lem) <- factors[,1], "must have the same number of columns")
    expect_error(featureLoadings(lem) <- loadings[,1], "must have the same number of columns")
    expect_error(factorData(lem) <- DataFrame(YAY=2), "one row per factor")
})

library(Matrix)
test_that("as.matrix works as expected", {
    expect_equal(as.matrix(lem), factors)

    alt <- rsparsematrix(nrow(lem), ncol(lem), density=0.1)
    sampleFactors(lem) <- alt
    expect_identical(sampleFactors(lem), alt)
    expect_equal(as.matrix(lem), as.matrix(alt))
})

test_that("manipulation of metadata is correct", {
    metadata(lem)$yay <- 1
    expect_equal(metadata(lem)$yay, 1)
    metadata(lem)$yay <- "stuff"
    expect_identical(metadata(lem)$yay, "stuff")
})
