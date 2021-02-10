# This tests the reduced.dim.matrix class
# library(testthat); library(SingleCellExperiment); source("setup.R"); source("test-reduced.dim.matrix.R")

pc <- matrix(runif(500), ncol=5)
rot <- matrix(rnorm(20), ncol=5)
rd.pc <- reduced.dim.matrix(pc, sdev=1:100, rotation=rot)

test_that("reduced.dim.matrix construction works as expected", {
    expect_identical(attr(rd.pc, "sdev"), 1:100)
    expect_identical(attr(rd.pc, "rotation"), rot)
    expect_identical(dim(rd.pc), dim(pc))

    stuff <- paste0("GENE_", seq_len(nrow(pc)))
    rownames(rd.pc) <- stuff
    expect_identical(rownames(rd.pc), stuff)

    expect_true(is.matrix(rd.pc))
    expect_true(is(rd.pc, "matrix"))
})

test_that("reduced.dim.matrix subsetting retains attributes", {
    sub.pc <- rd.pc[1:10,]
    expect_identical(attr(sub.pc, "sdev"), 1:100)
    expect_identical(attr(sub.pc, "rotation"), rot)
    attributes(sub.pc) <- attributes(sub.pc)[c('dim', 'dimnames')]
    expect_identical(sub.pc, pc[1:10,])

    sub.pc <- rd.pc[,1:2]
    expect_identical(attr(rd.pc, "sdev"), 1:100)
    expect_identical(attr(rd.pc, "rotation"), rot)
    attributes(sub.pc) <- attributes(sub.pc)[c('dim', 'dimnames')]
    expect_identical(sub.pc, pc[,1:2])

    sub.pc <- rd.pc[,1]
    expect_null(attr(sub.pc, "sdev"))
    expect_null(attr(sub.pc, "rotation"))
})

test_that("reduced.dim.matrix combining retains attributes", {
    com.pc <- rbind(rd.pc, rd.pc[1:10,])
    expect_identical(attr(com.pc, "sdev"), 1:100)
    expect_identical(attr(com.pc, "rotation"), rot)
    attributes(com.pc) <- attributes(com.pc)[c('dim', 'dimnames')]
    expect_identical(com.pc, rbind(pc, pc[1:10,]))

    com.pc <- cbind(rd.pc, rd.pc[,1:2])
    expect_identical(attr(com.pc, "sdev"), 1:100)
    expect_identical(attr(com.pc, "rotation"), rot)
    attributes(com.pc) <- attributes(com.pc)[c('dim', 'dimnames')]
    expect_identical(com.pc, cbind(pc, pc[,1:2]))

    # What happens with differences in attributes?
    expect_warning(com.pc <- rbind(rd.pc, pc[1:10,]), "mismatched") 
    expect_identical(com.pc, rbind(pc, pc[1:10,]))

    expect_warning(com.pc <- cbind(rd.pc, pc[,1:2]), "mismatched") 
    expect_identical(com.pc, cbind(pc, pc[,1:2]))
})

test_that("reduced.dim.matrix works in its intended environment", {
    sce <- loaded
    reducedDim(sce, "PCA") <- reduced.dim.matrix(reducedDim(sce, "PCA"), sdev=1:100)

    subsce <- sce[1:10,]
    expect_identical(attr(reducedDim(subsce, "PCA"), "sdev"), 1:100)

    comsce <- cbind(sce, sce)
    expect_identical(attr(reducedDim(comsce, "PCA"), "sdev"), 1:100)
})
