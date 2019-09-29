# Tests for proper functioning of the altExp methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-altExp.R")

sce <- loaded

test_that("altExp getters work correctly", {
    expect_identical(altExp(sce), se1)
    expect_identical(altExp(sce, 2), se2)
    expect_identical(altExp(sce, "Protein"), se2)

    expect_identical(altExps(sce, withColData=FALSE), List(Spike=se1, Protein=se2))
    expect_identical(altExpNames(sce), c("Spike", "Protein"))

    # Carries over colData.
    colData(sce)$stuff <- runif(ncol(sce))
    expect_identical(sce$stuff, altExp(sce, withColData=TRUE)$stuff)
    expect_identical(sce$stuff, altExp(sce, "Protein", withColData=TRUE)$stuff)
    expect_identical(NULL, altExp(sce, withColData=FALSE)$stuff)

    expect_identical(sce$stuff, altExps(sce, withColData=TRUE)[[1]]$stuff)
    expect_identical(NULL, altExps(sce, withColData=FALSE)[[1]]$stuff)
})

test_that("altExp setters work correctly", {
    se3 <- se2
    assay(se3) <- assay(se3) + 1
    altExp(sce) <- se3

    expect_identical(altExp(sce), se3)
    expect_identical(altExp(sce, 1), se3)
    expect_identical(altExp(sce, "Spike"), se3)

    altExp(sce, 2) <- se3
    expect_identical(altExp(sce, 2), se3)
    expect_identical(altExp(sce, "Protein"), se3)

    altExp(sce, 1) <- NULL
    expect_identical(altExpNames(sce), "Protein")
    altExp(sce, 1) <- NULL
    expect_identical(altExpNames(sce), character(0))

    altExp(sce) <- se2
    expect_identical(altExpNames(sce), "unnamed1")
    expect_identical(altExp(sce), se2)
    altExp(sce) <- se3
    expect_identical(altExpNames(sce), "unnamed1")
    expect_identical(altExp(sce), se3)

    expect_error(altExp(sce, 5) <- se3, "out of bounds")
})

test_that("altExps setters work correctly", {
    altExps(sce) <- NULL
    expect_identical(unname(altExps(sce)), List())
    altExps(sce) <- list(whee=se1, blah=se2)
    expect_identical(altExpNames(sce), c("whee", "blah"))
    expect_identical(altExp(sce,1), se1)

    # Works without names.
    altExps(sce) <- list(se1, se2)
    expect_identical(altExpNames(sce), c("unnamed1", "unnamed2"))
    expect_identical(altExp(sce,1), se1)
    expect_identical(altExp(sce,2), se2)
})

test_that("altExpNames setters work correctly", {
    altExpNames(sce) <- c("A", "B")
    expect_identical(altExpNames(sce), c("A", "B"))

    expect_error(altExpNames(empty) <- c("A", "B"), "more column names")
    expect_error(altExpNames(empty) <- NULL, "unable to find an inherited method")
    expect_error(altExpNames(sce) <- LETTERS, "more column names")
})

test_that("splitAltExps works correctly", {
    feat.type <- sample(c("endog", "ERCC", "CITE"), nrow(empty),
        replace=TRUE, p=c(0.8, 0.1, 0.1))
    out <- splitAltExps(empty, feat.type)

    expect_identical(assay(out), assay(empty[feat.type=="endog",]))
    expect_identical(altExp(out, "ERCC"), empty[feat.type=="ERCC",])
    expect_identical(altExp(out, "CITE"), empty[feat.type=="CITE",])

    # Handles alternative reference.
    out <- splitAltExps(empty, feat.type, ref="ERCC")

    expect_identical(assay(out), assay(empty[feat.type=="ERCC",]))
    expect_identical(altExp(out, "endog"), empty[feat.type=="endog",])
    expect_identical(altExp(out, "CITE"), empty[feat.type=="CITE",])

    # Clears out colData
    empty$blah <- sample(LETTERS, ncol(empty), replace=TRUE)
    out <- splitAltExps(empty, feat.type)
    expect_identical(colnames(colData(altExp(out, withColData=FALSE))), character(0))
})

test_that("swapAltExp works correctly", {
    feat.type <- sample(c("endog", "ERCC", "CITE"), nrow(empty),
        replace=TRUE, p=c(0.8, 0.1, 0.1))
    ref <- splitAltExps(empty, feat.type)
    ref$A <- seq_len(ncol(ref))

    swapped <- swapAltExp(ref, "CITE", save="all")
    expect_identical(assay(swapped), assay(altExp(ref, "CITE")))
    expect_identical(colData(swapped), colData(ref))
    expect_identical(altExpNames(swapped), c(altExpNames(ref), "all"))

    swapped2 <- swapAltExp(swapped, "all")
    expect_identical(assay(swapped2), assay(ref))
    expect_identical(colData(swapped2), colData(ref))
    expect_identical(altExps(swapped2), altExps(swapped))

    swapped3 <- swapAltExp(swapped, "CITE", withColData=FALSE)
    expect_identical(ncol(colData(swapped3)), 0L)
})

test_that("getters and setters throw appropriate errors", {
    expect_error(altExp(empty), "no available entries")
    expect_error(altExp(sce, 3), "subscript contains out-of-bounds indices")
    expect_error(altExp(sce, "dummy"), "invalid subscript")
    expect_error(altExp(sce, 3) <- se1, "out of bounds")
})

test_that(".precheck_altExp throws appropriate errors", {
    expect_error(
        SingleCellExperiment:::.precheck_altExp(sce, assay(sce)),
        "should be a SummarizedExperiment object"
    )

    expect_error(
        SingleCellExperiment:::.precheck_altExp(sce, SummarizedExperiment()),
        "should have the same number of columns as 'x'"
    )
})
