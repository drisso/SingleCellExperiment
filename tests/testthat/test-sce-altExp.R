# Tests for proper functioning of the altExp methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-altExp.R")

sce <- loaded

test_that("altExp getters work correctly", {
    expect_identical(altExp(sce), se1)
    expect_identical(altExp(sce, 2), se2)
    expect_identical(altExp(sce, "Protein"), se2)

    expect_identical(altExps(sce), List(Spike=se1, Protein=se2))
    expect_identical(altExpNames(sce), c("Spike", "Protein"))

    # Carries over column names.
    copy <- sce
    colnames(copy) <- seq_len(ncol(copy))
    expect_identical(colnames(altExp(copy)), colnames(copy))
    expect_null(colnames(altExp(copy, withDimnames=FALSE)))

    # Carries over colData.
    expect_identical(sce$sizeFactor, altExp(sce, withColData=TRUE)$sizeFactor)
    expect_identical(sce$sizeFactor, altExp(sce, "Protein", withColData=TRUE)$sizeFactor)
    expect_identical(NULL, altExp(sce, withColData=FALSE)$sizeFactor)

    expect_identical(sce$sizeFactor, altExps(sce, withColData=TRUE)[[1]]$sizeFactor)
    expect_identical(NULL, altExps(sce, withColData=FALSE)[[1]]$sizeFactor)

    # Prepends colData properly.
    copy <- sce
    u <- runif(ncol(sce))
    altExp(copy)$BLAH <- u
    expect_identical(colData(altExp(copy, withColData=TRUE)), cbind(colData(sce), colData(altExp(copy))))
    expect_identical(colData(altExp(copy)), DataFrame(BLAH=u))
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

test_that("altExp setter handles name mismatches", {
    colnames(se2) <- seq_len(ncol(se2))
    expect_warning(altExp(sce, 2) <- se2, "not the same")
    expect_warning(altExp(sce, 2, withDimnames=FALSE) <- se2, NA)

    colnames(sce) <- colnames(se2) <- seq_len(ncol(se2))
    expect_warning(altExp(sce, 2) <- se2, NA)
})

test_that("altExp setter correctly removes cloned colData.", {
    copy <- altExp(sce, 2, withColData=TRUE)

    # Attempt 1:
    expect_warning(altExp(sce, 2) <- copy, NA)
    expect_false(is.null(altExp(sce, 2)$sizeFactor))

    # Attempt 2:
    expect_warning(altExp(sce, 2, withColData=TRUE) <- copy, NA)
    expect_null(altExp(sce, 2)$sizeFactor)

    # Attempt 3:
    expect_warning(altExp(sce, 2, withColData=TRUE) <- se2, "left-most")
    expect_null(altExp(sce, 2)$sizeFactor)

    # Attempt 4:
    colData(copy) <- cbind(colData(copy), X=runif(ncol(copy)))
    expect_warning(altExp(sce, 2, withColData=TRUE) <- copy, NA)
    expect_null(altExp(sce, 2)$sizeFactor)
    expect_false(is.null(altExp(sce, 2)$X))

    # Attempt 5:
    expect_warning(rownames(altExp(sce, 2, withColData=TRUE)) <- 1:5, NA)
    expect_null(altExp(sce, 2)$sizeFactor)

    # Attempt 6:
    expect_warning(altExp(sce, 2, withColData=TRUE)$stuff <- 2, NA)
    expect_null(altExp(sce, 2)$sizeFactor)
    expect_true(all(altExp(sce, 2)$stuff==2))
})

test_that("altExps setters work correctly", {
    altExps(sce) <- NULL
    expect_identical(unname(altExps(sce)), List())
    altExps(sce) <- list(whee=se1, blah=se2)
    expect_identical(altExpNames(sce), c("whee", "blah"))
    expect_identical(altExp(sce,1), se1)

    # Works without names.
    expect_warning(altExps(sce) <- list(se1, se2), "NULL")
    expect_identical(altExpNames(sce), c("unnamed1", "unnamed2"))
    expect_identical(altExp(sce,1), se1)
    expect_identical(altExp(sce,2), se2)

    expect_warning(altExps(sce) <- list(X=se1, se2), "empty")
    expect_identical(altExpNames(sce), c("X", "unnamed1"))

    # Handles non-syntactical names.
    altExps(sce) <- list(`first thing`=se1, `second thing`=se2)
    expect_identical(altExpNames(sce), c("first thing", "second thing"))
})

test_that("altExps setter responds to withDimnames=", {
    colnames(se2) <- seq_len(ncol(se2))
    expect_warning(altExps(sce) <- list(YAY=se2), "not the same")
    expect_warning(altExps(sce, withDimnames=FALSE) <- list(YAY=se2), NA)

    colnames(sce) <- colnames(se2) <- seq_len(ncol(se2))
    expect_warning(altExps(sce) <- list(YAY=se2), NA)
})

test_that("altExps setter responds to withColData=", {
    copy <- altExps(sce, withColData=TRUE)
    expect_warning(altExps(sce) <- copy, NA)
    expect_false(is.null(altExp(sce)$sizeFactor))

    expect_warning(altExps(sce, withColData=TRUE) <- copy, NA)
    expect_null(altExp(sce)$sizeFactor)

    colData(copy[[2]]) <- cbind(stuff=runif(ncol(sce)), colData(copy[[2]]))
    expect_warning(altExps(sce, withColData=TRUE) <- copy, "left-most")
    expect_false(is.null(altExp(sce, 2)$stuff))
})

test_that("altExps getters/setters preserve mcols and metadata", {
    stuff <- List(whee=se1, blah=se2)
    mcols(stuff)$A <- c("one", "two")
    metadata(stuff)$B <- "three"

    altExps(sce) <- stuff
    out <- altExps(sce)
    expect_identical(mcols(out), mcols(stuff))
    expect_identical(metadata(out), metadata(stuff))
})

test_that("altExpNames setters work correctly", {
    altExpNames(sce) <- c("A", "B")
    expect_identical(altExpNames(sce), c("A", "B"))

    expect_warning(altExpNames(sce) <- c("X", ""), "empty")
    expect_identical(altExpNames(sce), c("X", "unnamed1"))

    expect_error(altExpNames(empty) <- c("A", "B"), "more column names")
    expect_error(altExpNames(empty) <- NULL, "unable to find an inherited method")
    expect_error(altExpNames(sce) <- LETTERS, "more column names")
})

test_that("splitAltExps works correctly", {
    feat.type <- sample(c("endog", "ERCC", "CITE"), nrow(empty),
        replace=TRUE, p=c(0.8, 0.1, 0.1))
    out <- splitAltExps(empty, feat.type)

    expect_identical(assay(out), assay(empty[feat.type=="endog",]))

    test <- empty
    mainExpName(test) <- NULL
    expect_identical(altExp(out, "ERCC"), test[feat.type=="ERCC",])
    expect_identical(altExp(out, "CITE"), test[feat.type=="CITE",])

    expect_identical(mainExpName(out), "endog")
    expect_null(mainExpName(altExp(out)))

    # Handles alternative references.
    out <- splitAltExps(empty, feat.type, ref="ERCC")

    expect_identical(assay(out), assay(empty[feat.type=="ERCC",]))
    expect_identical(altExp(out, "endog"), test[feat.type=="endog",])
    expect_identical(altExp(out, "CITE"), test[feat.type=="CITE",])
    expect_identical(mainExpName(out), "ERCC")

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

    swapped <- swapAltExp(ref, "CITE")
    expect_identical(assay(swapped), assay(altExp(ref, "CITE")))
    expect_identical(colData(swapped), colData(ref))

    expect_identical(altExpNames(swapped), c("ERCC", "endog"))
    expect_null(mainExpName(altExp(swapped, "endog")))
    expect_identical(mainExpName(swapped), "CITE")

    # Operation is perfectly reversible.
    swapped2 <- swapAltExp(swapped, "endog")
    altExps(swapped2) <- rev(altExps(swapped2))
    expect_identical(swapped2, ref)

    # More predictable behavior with withColData=FALSE.
    swapped3 <- swapAltExp(ref, "CITE", withColData=FALSE)
    expect_identical(ncol(colData(swapped3)), 0L)
    expect_identical(colnames(colData(altExp(swapped3, "endog"))), "A")
})

test_that("getters and setters throw appropriate errors", {
    expect_error(altExp(empty), "no available entries")
    expect_error(altExp(sce, 3), "subscript contains out-of-bounds indices")
    expect_error(altExp(sce, "dummy"), "invalid subscript")
    expect_error(altExp(sce, 3) <- se1, "out of bounds")
})

test_that(".precheck_altExp throws appropriate errors", {
    expect_error(
        altExp(sce) <- assay(sce),
        "SummarizedExperiment"
    )

    expect_error(
        altExp(sce) <- SummarizedExperiment(),
        "number of columns"
    )
})
