# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-combine.R")

sce <- loaded

test_that("rbind works correctly in the basic case", {
    shuffled <- sample(nrow(v))
    sce.alt <- sce[shuffled,]

    sce2 <- rbind(sce, sce.alt)
    expect_equivalent(assay(sce2), rbind(assay(sce), assay(sce.alt)))
    expect_identical(sizeFactors(sce2), sizeFactors(sce))
    expect_identical(isSpike(sce2), c(isSpike(sce), isSpike(sce.alt)))
    expect_identical(reducedDims(sce2), reducedDims(sce))
})

test_that("rbind respects the internal fields correctly", {
    # Respects the internal colData.
    sce2 <- sce
    int_colData(sce2)$X <- runif(ncol(sce2))
    sce3 <- rbind(sce, sce2)
    expect_identical(int_colData(sce3)$X, int_colData(sce2)$X)

    # Respects reordered internal elementMetadata
    int_elementMetadata(sce) <- cbind(int_elementMetadata(sce), DataFrame(A=runif(nrow(sce)), B=runif(nrow(sce))))
    alpha <- rbind(sce, sce)
    alt.sce <- sce
    int_elementMetadata(alt.sce) <- int_elementMetadata(alt.sce)[,ncol(int_elementMetadata(alt.sce)):1]
    bravo <- rbind(sce, alt.sce)
    expect_identical(alpha, bravo)
})

test_that("rbind handles errors in internal fields correctly", {
    sce2 <- sce
    int_colData(sce)$X <- runif(ncol(sce))
    int_colData(sce2)$X <- runif(ncol(sce2))
    expect_error(rbind(sce, sce2), "'int_colData'")

    # Throws errors upon mismatch in the internal elementMetadata.
    sce.err <- sce
    int_elementMetadata(sce.err)$X <- "YAY"
    expect_error(rbind(sce.err, sce), "'int_elementMetadata'")
})

test_that("cbind works correctly in the basic case", {
    shuffled <- sample(ncells)
    sce.alt <- sce[,shuffled]

    sce2 <- cbind(sce, sce.alt)
    expect_equivalent(assay(sce2), cbind(assay(sce), assay(sce.alt)))
    expect_identical(sizeFactors(sce2), c(sizeFactors(sce), sizeFactors(sce.alt)))
    expect_identical(reducedDim(sce2, "PCA"), rbind(reducedDim(sce, "PCA"), reducedDim(sce.alt, "PCA")))
    expect_identical(isSpike(sce2), isSpike(sce))
    expect_identical(altExp(sce2), cbind(altExp(sce), altExp(sce.alt)))
})

test_that("cbind respects the internal fields correctly", {
    # Respects the internal elementMetadata.
    sce2 <- sce
    int_elementMetadata(sce2)$X <- runif(nrow(sce2))
    sce3 <- cbind(sce, sce2)
    expect_identical(int_elementMetadata(sce3)$X, int_elementMetadata(sce2)$X)

    # Respects reordered internal colData.
    alpha <- cbind(sce, sce)
    alt.sce <- sce
    int_colData(alt.sce) <- int_colData(alt.sce)[,ncol(int_colData(alt.sce)):1]
    bravo <- cbind(sce, alt.sce)
    expect_identical(alpha, bravo)
    sce.err <- sce
    sizeFactors(sce.err) <- NULL
    expect_error(cbind(sce.err, sce), "'int_colData'")

    sce.err <- sce
    reducedDim(sce.err, "PCA") <- NULL
    expect_error(cbind(sce.err, sce), "'int_colData'")

    sce.err <- sce
    altExp(sce.err, 1) <- NULL
    expect_error(cbind(sce.err, sce), "'int_colData'")
})

