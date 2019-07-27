# Tests for proper functioning of the altExp methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-altExp.R")

sce <- loaded

test_that("altExp getters work correctly", {
    expect_identical(altExp(sce), se1)
    expect_identical(altExp(sce, 2), se2)
    expect_identical(altExp(sce, "Protein"), se2)
    expect_error(altExp(empty), "no alternative experiments")

    expect_identical(altExps(sce, withColData=FALSE), List(Spike=se1, Protein=se2))
    expect_identical(altExpNames(sce), c("Spike", "Protein"))

    # Carries over colData.
    colData(sce)$stuff <- runif(ncol(sce))
    expect_identical(sce$stuff, altExp(sce)$stuff)
    expect_identical(NULL, altExp(sce, withColData=FALSE)$stuff)

    expect_identical(sce$stuff, altExps(sce)[[1]]$stuff)
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
    expect_identical(altExpNames(sce), character(2))
    expect_identical(altExp(sce,1), se1)
    expect_identical(altExp(sce,2), se2)
})

test_that("altExpNames setters work correctly", {
    altExpNames(sce) <- c("A", "B")
    expect_identical(altExpNames(sce), c("A", "B"))

    expect_error(altExpNames(empty) <- c("A", "B"), "no alternative experiments")
    expect_error(altExpNames(empty) <- NULL, NA) # coerces the character.
    expect_error(altExpNames(sce) <- LETTERS, "more column names")
})

test_that("splitSCEByAlt works correctly", {
     feat.type <- sample(c("endog", "ERCC", "CITE"), nrow(empty),
         replace=TRUE, p=c(0.8, 0.1, 0.1))
     out <- splitSCEByAlt(empty, feat.type)

     expect_identical(assay(out), assay(empty[feat.type=="endog",]))
     expect_identical(altExp(out, "ERCC"), empty[feat.type=="ERCC",])
     expect_identical(altExp(out, "CITE"), empty[feat.type=="CITE",])

     # Handles alternative reference. 
     out <- splitSCEByAlt(empty, feat.type, ref="ERCC")

     expect_identical(assay(out), assay(empty[feat.type=="ERCC",]))
     expect_identical(altExp(out, "endog"), empty[feat.type=="endog",])
     expect_identical(altExp(out, "CITE"), empty[feat.type=="CITE",])
})
