# Checks for proper functioning of the methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-methods.R")

sce <- empty

test_that("size factor getters/setters are functioning", {
    sf1 <- 2^rnorm(ncells)
    sizeFactors(sce) <- sf1
    expect_identical(sizeFactors(sce), sf1)

    # Manual deletion.
    sizeFactors(sce) <- NULL
    expect_identical(sizeFactors(sce), NULL)

    # Additional actions work.
    expect_warning(sizeFactors(sce, onAbsence="warn"), "NULL")
    expect_error(sizeFactors(sce, onAbsence="error"), "NULL")
})

test_that("column label getters/setters are functioning", {
    labels <- sample(letters, ncol(sce), replace=TRUE) 
    colLabels(sce) <- labels
    expect_identical(colLabels(sce), labels)

    # Manual deletion.
    colLabels(sce) <- NULL
    expect_identical(colLabels(sce), NULL)

    # Additional actions work.
    expect_warning(colLabels(sce, onAbsence="warn"), "NULL")
    expect_error(colLabels(sce, onAbsence="error"), "NULL")
})

test_that("row subset getters/setters are functioning", {
    of.interest <- 1:10
    rowSubset(sce) <- of.interest
    expect_identical(which(rowSubset(sce)), of.interest)

    keep <- rbinom(nrow(sce), 1, 0.5)==1
    rowSubset(sce) <- keep
    expect_identical(rowSubset(sce), keep)

    rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
    labels <- sample(rownames(sce), ncol(sce), replace=TRUE) 
    rowSubset(sce) <- labels
    expect_identical(rowSubset(sce), rownames(sce) %in% labels)

    # Responds to non-default locations.
    of.interest <- 1:10
    rowSubset(sce, "hvgs") <- of.interest
    expect_identical(which(rowSubset(sce, "hvgs")), of.interest)

    # Manual deletion.
    rowSubset(sce) <- NULL
    expect_identical(rowSubset(sce), NULL)

    # Additional actions work.
    expect_warning(rowSubset(sce, onAbsence="warn"), "NULL")
    expect_error(rowSubset(sce, onAbsence="error"), "NULL")
})

test_that("object version extraction works", {
    expect_identical(objectVersion(sce), packageVersion("SingleCellExperiment"))
})

test_that("special colData/rowData getters/setters work", {
    int_elementMetadata(sce) <- DataFrame(STUFF=rbinom(nrow(v), 1, 0.2)==1)
    int_colData(sce) <- DataFrame(WHEE=2^rnorm(ncells))

    random_coldata <- DataFrame(a=rnorm(ncells), b=runif(ncells, 0, 1))
    colData(sce) <- random_coldata
    expect_identical(colData(sce, use.names=FALSE), random_coldata)
    expect_identical(colData(sce), colData(sce, internal=FALSE))
    expect_identical(colData(sce, internal=TRUE), cbind(colData(sce), int_colData(sce)))

    random_rowdata <- DataFrame(a=rnorm(NROW(sce)), b=runif(NROW(sce), 0, 1))
    rowData(sce) <- random_rowdata
    expect_identical(rowData(sce, use.names=FALSE), random_rowdata)
    expect_identical(rowData(sce), rowData(sce, internal=FALSE))
    expect_identical(rowData(sce, internal=TRUE), cbind(rowData(sce), int_elementMetadata(sce)))

    # Passes arguments correctly down the line.
    rout <- rowData(sce, use.names=FALSE)
    expect_identical(rownames(rout), NULL)
    cout <- colData(sce, use.names=FALSE)
    expect_identical(rownames(cout), NULL)

    sceN <- sce
    colnames(sceN) <- paste("Cell", seq_len(ncol(sceN)))
    rownames(sceN) <- paste("Cell", seq_len(nrow(sceN)))

    rout <- rowData(sceN, use.names=TRUE)
    expect_identical(rownames(rout), rownames(sceN))
    cout <- colData(sceN, use.names=TRUE)
    expect_identical(rownames(cout), colnames(sceN))

    # Warnings upon overlaps.
    rowData(sce)$STUFF <- rnorm(NROW(sce))
    expect_warning(rowData(sce, internal=TRUE), "overlapping names in internal and external rowData")

    colData(sce)$WHEE <- rnorm(ncells)
    expect_warning(colData(sce, internal=TRUE), "overlapping names in internal and external colData")
})

test_that("assay getters/setters work", {
    v2 <- matrix(runif(20000), ncol=ncells)
    counts(sce) <- v2
    expect_equivalent(counts(sce), v2)

    v3 <- log2(v2)
    logcounts(sce) <- v3
    expect_equivalent(counts(sce), v2)
    expect_equivalent(logcounts(sce), v3)

    cpm(sce) <- v3 + v2
    expect_equivalent(cpm(sce), v3+v2)
    tpm(sce) <- v3 - v2
    expect_equivalent(tpm(sce), v3-v2)

    v4 <- v2 * v3
    weights(sce) <- v4
    expect_equivalent(weights(sce), v4)

    counts(sce) <- NULL
    expect_equivalent(logcounts(sce), v3)
    expect_error(counts(sce), "invalid subscript")
})

test_that("assay getters/setters respect withDimnames", {
    sce_rownames <- paste0("G", seq_len(20000 / ncells))
    sce_colnames <- paste0("S", seq_len(ncells))
    v2 <- matrix(runif(20000), ncol=ncells)
    counts(sce) <- v2
    rownames(sce) <- sce_rownames
    colnames(sce) <- sce_colnames
    expect_identical(dimnames(counts(sce)), list(sce_rownames, sce_colnames))
    expect_identical(dimnames(counts(sce, withDimnames=FALSE)), NULL)

    v3 <- log2(v2)
    logcounts(sce, withDimnames=FALSE) <- v3
    expect_identical(dimnames(logcounts(sce)), list(sce_rownames, sce_colnames))
    expect_identical(dimnames(logcounts(sce, withDimnames=FALSE)), NULL)

    cpm(sce, withDimnames=FALSE) <- v3 + v2
    expect_identical(dimnames(cpm(sce)), list(sce_rownames, sce_colnames))
    expect_identical(dimnames(cpm(sce, withDimnames=FALSE)), NULL)

    tpm(sce, withDimnames=FALSE) <- v3 - v2
    expect_equivalent(tpm(sce), v3-v2)
    expect_identical(dimnames(tpm(sce)), list(sce_rownames, sce_colnames))
    expect_identical(dimnames(tpm(sce, withDimnames=FALSE)), NULL)
})
