# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-combine.R")

sce <- loaded

test_that("rbind works correctly in the basic case", {
    shuffled <- sample(nrow(v))
    sce.alt <- sce[shuffled,]

    sce2 <- rbind(sce, sce.alt)
    expect_equivalent(assay(sce2), rbind(assay(sce), assay(sce.alt)))
    expect_identical(sizeFactors(sce2), sizeFactors(sce))
    expect_identical(reducedDims(sce2), reducedDims(sce))

    dual1 <- SingleCellExperiment:::DualSubset(rowPair(sce))
    dual2 <- SingleCellExperiment:::DualSubset(rowPair(sce.alt))
    expect_identical(
        SingleCellExperiment:::DualSubset(rowPair(sce2)),
        c(dual1, dual2) 
    )
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

    # Don't concatenate names when merging metadata().
    sce4 <- rbind(A=sce, B=sce)
    expect_identical(objectVersion(sce4), objectVersion(sce))
})

test_that("cbind works correctly in the basic case", {
    shuffled <- sample(ncells)
    sce.alt <- sce[,shuffled]

    sce2 <- cbind(sce, sce.alt)
    expect_equivalent(assay(sce2), cbind(assay(sce), assay(sce.alt)))
    expect_identical(sizeFactors(sce2), c(sizeFactors(sce), sizeFactors(sce.alt)))
    expect_identical(reducedDim(sce2, "PCA"), rbind(reducedDim(sce, "PCA"), reducedDim(sce.alt, "PCA")))
    expect_identical(altExp(sce2), cbind(altExp(sce), altExp(sce.alt)))

    dual1 <- SingleCellExperiment:::DualSubset(colPair(sce))
    dual2 <- SingleCellExperiment:::DualSubset(colPair(sce.alt))
    expect_identical(
        SingleCellExperiment:::DualSubset(colPair(sce2)),
        c(dual1, dual2) 
    )
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

    alt.sce <- sce
    reducedDims(alt.sce) <- rev(reducedDims(alt.sce))
    altExps(alt.sce) <- rev(altExps(alt.sce))
    bravo <- cbind(sce, alt.sce)
    expect_identical(alpha, bravo)
})

test_that("cbind handles errors in the internal fields correctly", {
    # Chokes correctly when presented with errors.
    sce.err <- sce
    reducedDim(sce.err, "PCA") <- NULL
    expect_error(cbind(sce.err, sce), "'int_colData'")

    sce.err <- sce
    altExp(sce.err, 1) <- NULL
    expect_error(cbind(sce.err, sce), "'int_colData'")

    # Don't concatenate names when merging metadata().
    sce4 <- rbind(A=sce, B=sce) 
    expect_identical(objectVersion(sce4), objectVersion(sce))
})

test_that("combineCols passes on delay when combining altExps", {
    sce.combined <- combineCols(sce, sce, use.names=FALSE, delayed=FALSE)
    expect_is(assay(altExp(sce.combined)), "matrix")

    sce.combined2 <- combineCols(sce, sce, use.names=FALSE, delayed=TRUE)
    expect_is(assay(altExp(sce.combined2)), "DelayedMatrix")
})

test_that("combineCols passes on fill when combining altExps", {
    sce1 <- sce[1:6, 1:10]
    rownames(sce1) <- LETTERS[seq_len(nrow(sce1))]
    colnames(sce1) <- letters[seq_len(ncol(sce1))]
    altExp(sce1) <- altExp(sce1)[1:6, ]

    sce2 <- sce[1:14, 11:20]
    rownames(sce2) <- LETTERS[1:14]
    colnames(sce2) <- letters[seq_len(ncol(sce2)) + ncol(sce1)]
    altExp(sce2) <- altExp(sce2)[4:7, ]

    sce.combined1 <- combineCols(sce1, sce2, use.names=TRUE)
    expect_true(anyNA(assay(sce.combined1)))

    sce.combined2 <- combineCols(sce1, sce2, use.names=TRUE, fill = -1)
    expect_true(any(assay(sce.combined2) == -1))

    sce.combined3 <- combineCols(sce1, sce2, use.names=TRUE)
    expect_true(anyNA(assay(altExp(sce.combined3))))

    sce.combined4 <- combineCols(sce1, sce2, use.names=TRUE, fill = -1)
    expect_true(any(assay(sce.combined4) == -1))
})
