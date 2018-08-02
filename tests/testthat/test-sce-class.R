# Checks for proper construction and get/setting of the slots of a SingleCellExperiment.
# library(SingleCellExperiment); library(testthat); source("test-sce-class.R")
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

    expect_equal(names(sce@reducedDims), character(0))
})

test_that("coercion from other classes works correctly", {
    sce <- SingleCellExperiment(u, rowData=rd, colData=cd)

    # Coercion from SummarizedExperiment
    se <- SummarizedExperiment(u, rowData = rd, colData = cd)
    sce2 <- as(se, "SingleCellExperiment")
    expect_equal(sce, sce2)

    # Checking that reduced dim names are set correctly.
    expect_equal(names(sce2@reducedDims), character(0))

    # Coercion from RangedSummarizedExperiment
    ranges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                      strand=sample(c("+", "-"), 200, TRUE))
    mcols(ranges) <- rd

    rse <- SummarizedExperiment(u, colData = cd, rowRanges = ranges)
    sce3 <- as(rse, "SingleCellExperiment")
    expect_equal(assays(sce), assays(sce3))
    expect_equal(rowData(sce), rowData(sce3))
    expect_equal(colData(sce), colData(sce3))
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
    expect_identical(nrow(SingleCellExperiment:::int_elementMetadata(sce)), nrow(sce))
    expect_identical(nrow(SingleCellExperiment:::int_colData(sce)), ncol(sce))
    expect_identical(length(SingleCellExperiment:::int_metadata(sce)), 3L)

    rextra <- rnorm(nrow(v))
    SingleCellExperiment:::int_elementMetadata(sce)$whee <- rextra
    expect_equal(rextra, SingleCellExperiment:::int_elementMetadata(sce)$whee)
    SingleCellExperiment:::int_elementMetadata(sce) <- DataFrame(1:5)
    expect_error(validObject(sce), "'nrow' of internal 'rowData' not equal to 'nrow(object)'", fixed=TRUE)

    cextra <- rnorm(ncells)
    SingleCellExperiment:::int_colData(sce)$stuff <- cextra
    expect_equal(cextra, SingleCellExperiment:::int_colData(sce)$stuff)
    SingleCellExperiment:::int_colData(sce) <- DataFrame(1:5)
    expect_error(validObject(sce), "'nrow' of internal 'colData' not equal to 'ncol(object)'", fixed=TRUE)

    SingleCellExperiment:::int_metadata(sce)$urg <- "I was here"
    expect_identical(SingleCellExperiment:::int_metadata(sce)$urg, "I was here")
})
