# Checks for proper construction and get/setting of the slots of a SingleCellExperiment.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-class.R")
context("SingleCellExperiment class")

set.seed(1000)
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
u <- matrix(rpois(20000, 5), ncol=ncells)
w <- matrix(runif(20000), ncol=ncells)
rd <- DataFrame(stuff=runif(nrow(v)))
cd <- DataFrame(whee=runif(ncells))

test_that("construction of the SCE works correctly", {
    sce <- SingleCellExperiment(assay=v)
    expect_equivalent(assay(sce), v)

    sce <- SingleCellExperiment(assay=list(counts=u, exprs=v))
    expect_equivalent(assay(sce, "counts"), u)
    expect_equivalent(assay(sce, "exprs"), v)

    assay(sce, "exprs") <- w
    expect_equivalent(assay(sce, "exprs"), w)
})

test_that("coercion from other classes works correctly", {
    # Coercion from SummarizedExperiment
    se <- SummarizedExperiment(u, rowData = rd, colData = cd)
    sce2 <- as(se, "SingleCellExperiment")
    expect_equal(sce2, SingleCellExperiment(u, rowData=rd, colData=cd)) # equality not identity due to environment in 'assays'.
    expect_true(validObject(sce2))

    # Coercion from RangedSummarizedExperiment
    ranges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                      strand=sample(c("+", "-"), 200, TRUE))
    mcols(ranges) <- rd

    rse <- SummarizedExperiment(u, colData = cd, rowRanges = ranges)
    sce3 <- as(rse, "SingleCellExperiment")
    expect_equal(sce3, SingleCellExperiment(u, colData=cd, rowRanges=ranges))
    expect_true(validObject(sce3))
})

test_that("manipulation of SE metadata is correct", {
    sce <- SingleCellExperiment(u, rowData=rd, colData=cd)

    # Checking colData and rowData
    expect_equal(rd, rowData(sce))
    expect_equal(cd, colData(sce))

    # Adding extra metadata fields.
    cextra <- rnorm(ncells)
    sce$blah <- cextra
    expect_equal(cextra, colData(sce)$blah)
    rextra <- rnorm(nrow(v))
    rowData(sce)$blah <- rextra
    expect_equal(rextra, rowData(sce)$blah)

    # Adding metadata fields.
    sce <- SingleCellExperiment(u, metadata=list(yay=1))
    expect_equal(metadata(sce)$yay, 1)
    metadata(sce)$yay <- "stuff"
    expect_identical(metadata(sce)$yay, "stuff")
})

test_that("internal functions work correctly", {
    sce <- SingleCellExperiment(assay=u)
    expect_identical(nrow(int_elementMetadata(sce)), nrow(sce))
    expect_identical(nrow(int_colData(sce)), ncol(sce))
    expect_identical(length(int_metadata(sce)), 1L)

    rextra <- rnorm(nrow(v))
    int_elementMetadata(sce)$whee <- rextra
    expect_equal(rextra, int_elementMetadata(sce)$whee)
    int_elementMetadata(sce) <- int_elementMetadata(sce)[1:5,]
    expect_error(validObject(sce), "'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'", fixed=TRUE)

    cextra <- rnorm(ncells)
    int_colData(sce)$stuff <- cextra
    expect_equal(cextra, int_colData(sce)$stuff)
    int_colData(sce) <- DataFrame(1:5)
    expect_error(validObject(sce), "'nrow' of 'int_colData' not equal to 'ncol(object)'", fixed=TRUE)

    int_metadata(sce)$urg <- "I was here"
    expect_identical(int_metadata(sce)$urg, "I was here")
})

test_that(".sce_show works", {
    expect_null(show(loaded))
})

