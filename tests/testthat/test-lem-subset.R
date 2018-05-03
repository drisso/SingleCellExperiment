# Checks the subsetting methods.
# library(SingleCellExperiment); library(testthat); source("test-lem-subset.R")
context("LEM subsetting")

set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
fdata <- DataFrame(WHEE=sample(LETTERS, ncol(factors)))
lem <- LinearEmbeddingMatrix(factors, loadings, fdata)

test_that("subsetting works correctly for different index types", {
    keep.dimnames <- FALSE
    for (x in 1:3) {
        if (x==1) {
            by.row <- sample(nrow(factors), nrow(factors)/2)
            by.col <- sample(ncol(factors), ncol(factors)/2)
        } else if (x==2) {
            by.row <- rbinom(nrow(factors), 1, 0.5)==1
            by.col <- rbinom(ncol(factors), 1, 0.5)==1
        } else {
            colnames(lem) <- paste0("Factor_", seq_len(ncol(factors)))
            rownames(lem) <- paste0("Cell_", seq_len(nrow(factors)))
            dimnames(factors) <- dimnames(lem)
            colnames(loadings) <- rownames(fdata) <- colnames(lem)

            by.row <- sample(rownames(lem), nrow(factors)/2)
            by.col <- sample(colnames(lem), ncol(factors)/2)
            keep.dimnames <- TRUE
        }

        # By row.
        lem.alt <- lem[by.row,]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), factors[by.row,])
        expect_identical(featureLoadings(lem.alt), loadings)
        expect_identical(factorData(lem.alt), fdata)

        # By column.
        lem.alt <- lem[,by.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), factors[,by.col])
        expect_identical(featureLoadings(lem.alt), loadings[,by.col])
        expect_identical(factorData(lem.alt), fdata[by.col,,drop=FALSE])

        # By row and column.
        lem.alt <- lem[by.row, by.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), factors[by.row,by.col])
        expect_identical(featureLoadings(lem.alt), loadings[,by.col])
        expect_identical(factorData(lem.alt), fdata[by.col,,drop=FALSE])
    }
})

test_that("subsetting works correctly with drop=TRUE", {
    # By row, with and without drop.
    keeper <- lem[1,]
    expect_identical(keeper, factors[1,])

    nodrop <- lem[1,,drop=FALSE]
    expect_identical(sampleFactors(nodrop, withDimnames=FALSE), factors[1,,drop=FALSE])
    expect_identical(featureLoadings(nodrop), loadings)
    expect_identical(factorData(nodrop), fdata)

    # By column, with and without drop.
    keeper <- lem[,1]
    expect_identical(keeper, factors[,1])

    nodrop <- lem[,1,drop=FALSE]
    expect_identical(sampleFactors(nodrop, withDimnames=FALSE), factors[,1,drop=FALSE])
    expect_identical(featureLoadings(nodrop), loadings[,1,drop=FALSE])
    expect_identical(factorData(nodrop), fdata[1,,drop=FALSE])

    # Handles names correctly.
    colnames(lem) <- paste0("Factor_", seq_len(ncol(factors)))
    rownames(lem) <- paste0("Cell_", seq_len(nrow(factors)))
    expect_identical(lem[,1], sampleFactors(lem)[,1])
    expect_identical(lem[2,], sampleFactors(lem)[2,])
    expect_identical(lem[1,2], sampleFactors(lem)[1,2])

    # Throws errors correctly.
    expect_error(lem[nrow(lem)+1,], "subscript out of bounds", fixed=TRUE)
    expect_error(lem["A",], "index out of bounds", fixed=TRUE)
})

test_that("subsetting assignment works correctly", {
    keep.dimnames <- FALSE
    for (x in 1:3) {
        if (x==1) {
            dest.row <- sample(nrow(factors), nrow(factors)/2)
            src.row <- sample(nrow(factors), nrow(factors)/2)
            dest.col <- sample(ncol(factors), ncol(factors)/2)
            src.col <- sample(ncol(factors), ncol(factors)/2)
        } else if (x==2) {
            dest.row <- rbinom(nrow(factors), 1, 0.5)==1
            src.row <- sample(dest.row)
            dest.col <- rbinom(ncol(factors), 1, 0.5)==1
            src.col <- sample(dest.col)
        } else {
            colnames(lem) <- paste0("Factor_", seq_len(ncol(factors)))
            rownames(lem) <- paste0("Gene_", seq_len(nrow(factors)))
            dimnames(factors) <- dimnames(lem)
            colnames(loadings) <- rownames(fdata) <- colnames(lem)

            dest.row <- sample(rownames(lem), nrow(factors)/2)
            src.row <- sample(rownames(lem), nrow(factors)/2)
            dest.col <- sample(colnames(lem), ncol(factors)/2)
            src.col <- sample(colnames(lem), ncol(factors)/2)
            keep.dimnames <- TRUE
        }

        # By row.
        lem.alt <- lem
        lem.alt[dest.row,] <- lem[src.row,]

        ref <- factors
        ref[dest.row,] <- ref[src.row,]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), ref)
        expect_identical(featureLoadings(lem.alt), loadings)
        expect_identical(factorData(lem.alt), fdata)

        # By column.
        lem.alt <- lem
        lem.alt[,dest.col] <- lem[,src.col]

        ref_sf <- factors
        ref_sf[,dest.col] <- ref_sf[,src.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), ref_sf)

        ref_fl <- loadings
        ref_fl[,dest.col] <- ref_fl[,src.col]
        expect_identical(featureLoadings(lem.alt), ref_fl)

        ref_fd <- fdata
        ref_fd[dest.col,] <- fdata[src.col,]
        expect_identical(factorData(lem.alt), ref_fd)

        # By row and column.
        lem.alt <- lem
        lem.alt[dest.row,dest.col] <- lem[src.row,src.col]

        ref_sf <- factors
        ref_sf[dest.row,dest.col] <- ref_sf[src.row,src.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), ref_sf)

        ref_fl <- loadings
        ref_fl[,dest.col] <- ref_fl[,src.col]
        expect_identical(featureLoadings(lem.alt), ref_fl)

        ref_fd <- fdata
        ref_fd[dest.col,] <- fdata[src.col,]
        expect_identical(factorData(lem.alt), ref_fd)
    }
})

test_that("subsetting works correctly with names", {
    rownames(lem) <- paste0("CELL", seq_len(nrow(lem)))
    colnames(lem) <- paste0("FACTOR", seq_len(ncol(lem)))
    out <- lem[,1:5]
    expect_identical(colnames(out), colnames(lem)[1:5])
    expect_identical(rownames(out), rownames(lem))

    out <- lem[1:5,]
    expect_identical(rownames(out), rownames(lem)[1:5])
    expect_identical(colnames(out), colnames(lem))

    # Repeating with subset assignment; names should NOT change.
    lem2 <- lem
    lem2[1:5,] <- lem[5:1,]
    expect_identical(rownames(lem2), rownames(lem))
    expect_identical(colnames(lem2), colnames(lem))

    lem2 <- lem
    lem2[,1:5] <- lem[,5:1]
    expect_identical(colnames(lem2), colnames(lem))
    expect_identical(rownames(lem2), rownames(lem))

    unnamed <- lem
    rownames(unnamed) <- NULL
    lem2 <- lem
    lem2[1:5,] <- unnamed[5:1,]
    expect_identical(rownames(lem2), rownames(lem))

    unnamed[1:5,] <- lem[5:1,]
    expect_identical(rownames(unnamed), NULL)
})


