# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-combine.R")

sce <- loaded

test_that("rbind works correctly", {
    shuffled <- sample(nrow(v))
    sce.alt <- sce[shuffled,]

    sce2 <- rbind(sce, sce.alt)
    expect_equivalent(assay(sce2), rbind(assay(sce), assay(sce.alt)))
    expect_identical(isSpike(sce2), c(isSpike(sce), isSpike(sce.alt)))
    expect_identical(sizeFactors(sce2), sizeFactors(sce))
    expect_identical(reducedDims(sce2), reducedDims(sce))

    sce.err <- sce
    isSpike(sce.err, "ERCC") <- NULL
    expect_error(rbind(sce.err, sce), "DataFrame 1 does not have 'is_spike_ERCC', 'is_spike'")
    expect_error(rbind(sce, sce.err), "DataFrame 2 does not have 'is_spike_ERCC', 'is_spike'")

    sce.lost <- sce
    sizeFactors(sce.lost) <- NULL
    sce2 <- rbind(sce.lost, sce)
    expect_identical(sizeFactors(sce2), NULL)
    sce2 <- rbind(sce, sce.lost)
    expect_identical(sizeFactors(sce2), sizeFactors(sce))

    sce.lost <- sce
    reducedDims(sce.lost) <- NULL
    sce2 <- rbind(sce.lost, sce)
    expect_equivalent(reducedDims(sce2), SimpleList())
    sce2 <- rbind(sce, sce.lost)
    expect_identical(reducedDims(sce2), reducedDims(sce))
})

test_that("cbind works correctly", {
    shuffled <- sample(ncells)
    sce.alt <- sce[,shuffled]

    sce2 <- cbind(sce, sce.alt)
    expect_equivalent(assay(sce2), cbind(assay(sce), assay(sce.alt)))
    expect_identical(sizeFactors(sce2), c(sizeFactors(sce), sizeFactors(sce.alt)))
    expect_identical(reducedDim(sce2, "PCA"), rbind(reducedDim(sce, "PCA"), reducedDim(sce.alt, "PCA")))
    expect_identical(isSpike(sce2), isSpike(sce))

    sce.err <- sce
    sizeFactors(sce.err) <- NULL
    expect_error(cbind(sce.err, sce), "DataFrame 1 does not have 'size_factor'")
    expect_error(cbind(sce, sce.err), "DataFrame 2 does not have 'size_factor'")

    sce.err <- sce
    reducedDim(sce.err, "PCA") <- NULL
    expect_error(cbind(sce.err, sce), "object 1 does not have 'PCA' in 'reducedDims'")
    expect_error(cbind(sce, sce.err), "object 2 does not have 'PCA' in 'reducedDims'")

    sce.err <- sce
    altExperiment(sce.err, 1) <- NULL
    expect_error(cbind(sce.err, sce), "object 1 does not have 'Spike' in 'altExperiments'")
    expect_error(cbind(sce, sce.err), "object 2 does not have 'Spike' in 'altExperiments'")

    sce.lost <- sce
    isSpike(sce.lost, "ERCC") <- NULL
    sce2 <- cbind(sce.lost, sce)
    expect_identical(isSpike(sce2), NULL)
    sce2 <- cbind(sce, sce.lost)
    expect_identical(isSpike(sce2), isSpike(sce))
})
