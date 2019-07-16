# Tests for proper functioning of the altExp methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-altExp.R")

sce <- loaded

test_that("altExp getters work correctly", {
    expect_identical(altExp(sce), se1)
    expect_identical(altExp(sce, 2), se2)
    expect_identical(altExp(sce, "Protein"), se2)
    expect_error(altExp(empty), "no alternative experiments")

    expect_identical(altExps(sce), List(Spike=se1, Protein=se2))
    expect_identical(altExpNames(sce), c("Spike", "Protein"))

    expect_identical(altAssay(sce), assay(se1))
    expect_identical(altAssay(sce, i=2), assay(se1, 2))
    expect_identical(altAssay(sce, "Protein"), assay(se2))
    expect_identical(altAssay(sce, 2, 2), assay(se2, 2))

    expect_identical(altAssayNames(sce), c("counts", "logcounts"))
    expect_identical(altAssayNames(sce, 2), c("counts", "logcounts"))
    expect_identical(altAssayNames(sce, "Protein"), c("counts", "logcounts"))
           
    expect_identical(altRowData(sce), rowData(se1))
    expect_identical(altRowData(sce, 2), rowData(se2))
    expect_identical(altRowData(sce, "Protein"), rowData(se2))

    expect_identical(altRowNames(sce), rownames(se1))
    expect_identical(altRowNames(sce, 2), rownames(se2))
    expect_identical(altRowNames(sce, "Protein"), rownames(se2))
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
})

test_that("altExpNames setters work correctly", {
    altExpNames(sce) <- c("A", "B")
    expect_identical(altExpNames(sce), c("A", "B"))

    expect_error(altExpNames(empty) <- c("A", "B"), "no alternative experiments")
    expect_error(altExpNames(sce) <- LETTERS, "more column names")
})

test_that("altAssay setters work correctly", {
    altAssay(sce) <- assay(se1) + 1
    expect_identical(altAssay(sce), assay(se1) + 1)
    altAssay(sce, 2) <- assay(se2) + 2
    expect_identical(altAssay(sce, 2), assay(se2) + 2)
    altAssay(sce, "Protein") <- assay(se2) + 3
    expect_identical(altAssay(sce, "Protein"), assay(se2) + 3)

    altAssay(sce, i=2) <- assay(se1, i=2) + 1
    expect_identical(altAssay(sce, i=2), assay(se1, 2) + 1)
    altAssay(sce, "Protein", i=2) <- assay(se2, i=2) + 1
    expect_identical(altAssay(sce, "Protein", i=2), assay(se2, 2) + 1)
})

test_that("altAssayNames setters work correctly", {
    altAssayNames(sce) <- c("Whee", "Boo")
    expect_identical(altAssayNames(sce), c("Whee", "Boo"))
    altAssayNames(sce, 2) <- c("Whee", "Boo")
    expect_identical(altAssayNames(sce, 2), c("Whee", "Boo"))
})

test_that("altRowData setters work correctly", {
    stuff <- paste0("X", altRowData(sce)$stuff)
    altRowData(sce)$stuff <- stuff
    expect_identical(altRowData(sce)$stuff, stuff)

    blah <- paste0("X", altRowData(sce, 2)$blah)
    altRowData(sce, 2)$blah <- blah
    expect_identical(altRowData(sce, 2)$blah, blah)

    blah2 <- paste0("Z", altRowData(sce, 2)$blah)
    altRowData(sce, "Protein")$blah <- blah2
    expect_identical(altRowData(sce, "Protein")$blah, blah2)
})

test_that("altRowNames setters work correctly", {
    rn <- altRowNames(sce) 
    rn <- paste0("X", rn)
    altRowNames(sce) <- rn
    expect_identical(altRowNames(sce), rn)

    rn <- altRowNames(sce, 2) 
    rn <- paste0("Y", rn)
    altRowNames(sce, 2) <- rn
    expect_identical(altRowNames(sce, 2), rn)

    rn <- altRowNames(sce, "Protein")
    rn <- paste0("Z", rn)
    altRowNames(sce, "Protein") <- rn
    expect_identical(altRowNames(sce, 2), rn)
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
