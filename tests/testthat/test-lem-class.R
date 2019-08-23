# Checks for proper construction and get/setting of the slots of a LinearEmbeddingMatrix.
# library(SingleCellExperiment); library(testthat); source("test-lem-class.R")
context("LinearEmbeddingMatrix class")
set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
lem <- LinearEmbeddingMatrix(factors, loadings)

test_that("LEM construction works correctly", {
    expect_equivalent(sampleFactors(lem, withDimnames=FALSE), factors)
    expect_identical(featureLoadings(lem), loadings)
    expect_identical(nrow(factorData(lem)), ncol(lem))
    expect_identical(ncol(factorData(lem)), 0L)

    # Trying again with an initialized factorData.
    fd <- DataFrame(YAY=seq_len(ncol(factors)))
    lem2 <- LinearEmbeddingMatrix(factors, loadings, factorData=fd)
    expect_identical(factorData(lem2), fd)
    expect_identical(lem2$YAY, fd$YAY)

    # Should throw errors if it doesn't make sense.
    expect_error(LinearEmbeddingMatrix(factors, loadings[,1:2]), "must have the same number of\n *columns")
    expect_error(LinearEmbeddingMatrix(factors, loadings, DataFrame(YAY=2)), "one row per factor")

    # Add metadata.
    md <- list(a = 1, b = "two")
    lem3 <- LinearEmbeddingMatrix(factors, loadings, factorData=fd, metadata=md)
    expect_identical(md, metadata(lem3))
})

test_that("LEM setters work correctly", {
    lem2 <- lem
    sampleFactors(lem2) <- factors * 2
    expect_identical(sampleFactors(lem2, withDimnames=FALSE), factors * 2)

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
    expect_error(sampleFactors(lem) <- factors[,1], "matrix-like")
    expect_error(featureLoadings(lem) <- loadings[,1], "matrix-like")
    expect_error(factorData(lem) <- 2, "DataFrame")
    expect_error(factorData(lem) <- DataFrame(YAY=2), "one row per factor")
})

test_that("getting with names makes sense", {
    rn <- paste0("Cell_", seq_len(nrow(lem)))
    rownames(lem) <- rn
    cn <- paste0("Factor_", seq_len(ncol(lem)))
    colnames(lem) <- cn

    expect_identical(rn, rownames(lem))
    expect_identical(cn, colnames(lem))

    expect_identical(sampleFactors(lem, withDimnames=FALSE), factors)
    ref <- factors
    dimnames(ref) <- list(rn, cn)
    expect_identical(sampleFactors(lem), ref)

    expect_identical(featureLoadings(lem, withDimnames=FALSE), loadings)
    ref <- loadings
    colnames(ref) <- cn
    expect_identical(featureLoadings(lem), ref)

    ref <- new("DFrame", nrows=ncol(lem))
    rownames(ref) <- cn
    expect_identical(factorData(lem, withDimnames=FALSE), ref)
    expect_identical(factorData(lem), ref) # Doesn't make a difference, as this is the reference location.
})

test_that("setting with names makes sense", {
    rn <- paste0("Cell_", seq_len(nrow(lem)))
    rownames(lem) <- rn
    cn <- paste0("Factor_", seq_len(ncol(lem)))
    colnames(lem) <- cn

    lem2 <- lem
    sampleFactors(lem2) <- factors * 2
    expect_identical(rownames(sampleFactors(lem2)), rn)
    expect_identical(colnames(sampleFactors(lem2)), cn)
    expect_equal(sampleFactors(lem2), sampleFactors(lem)*2)

    lem2 <- lem
    featureLoadings(lem2) <- loadings * 2
    expect_identical(colnames(featureLoadings(lem2)), cn)
    expect_equal(featureLoadings(lem2), featureLoadings(lem)*2)

    # This will delete the column names.
    lem2 <- lem
    factorData(lem2) <- new("DFrame", nrows=ncol(lem))
    expect_equal(ncol(factorData(lem2)), 0L)
    expect_identical(rownames(factorData(lem2)), NULL)

    # Throwing errors when the names don't match up.
    factor2 <- factors
    rownames(factor2) <- paste("X", rn)
    expect_error(sampleFactors(lem) <- factor2, "must match")
    factor2 <- factors
    colnames(factor2) <- paste("X", cn)
    expect_error(sampleFactors(lem) <- factor2, "must match")

    loading2 <- loadings
    colnames(loading2) <- paste0("Y", cn)
    expect_error(featureLoadings(lem) <- loading2, "must match")

    ref <- new("DFrame", nrows=ncol(lem))
    rownames(ref) <- paste0("Y", cn)
    expect_error(factorData(lem) <- ref, "must match")
})

library(Matrix)
test_that("as.matrix works as expected", {
    factors2 <- factors
    dimnames(factors2) <- dimnames(lem)
    expect_equal(as.matrix(lem), factors2)

    alt <- rsparsematrix(nrow(lem), ncol(lem), density=0.1)
    sampleFactors(lem) <- alt
    expect_identical(sampleFactors(lem, withDimnames=FALSE), alt)
    expect_equal(as.matrix(lem), as.matrix(alt))
})

test_that("manipulation of metadata is correct", {
    metadata(lem)$yay <- 1
    expect_equal(metadata(lem)$yay, 1)
    metadata(lem)$yay <- "stuff"
    expect_identical(metadata(lem)$yay, "stuff")
})

test_that("show method works", {
    expect_output(show(LinearEmbeddingMatrix()))
})
