# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("test-lem-subset.R")

set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
fdata <- DataFrame(WHEE=sample(LETTERS, ncol(factors)))
lem <- LinearEmbeddingMatrix(factors, loadings, fdata)

test_that("subsetting works correctly", {
    # By row.
    shuffled <- sample(nrow(lem))
    lem.alt <- lem[shuffled,]

    expect_identical(sampleFactors(lem.alt), factors[shuffled,])
    expect_identical(featureLoadings(lem.alt), loadings)
    expect_identical(factorData(lem.alt), fdata)

    # By column.
    shuffled <- sample(ncol(factors))
    lem.alt <- lem[,shuffled]

    expect_identical(sampleFactors(lem.alt), factors[,shuffled])
    expect_identical(featureLoadings(lem.alt), loadings[,shuffled])
    expect_identical(factorData(lem.alt), fdata[shuffled,,drop=FALSE])

    # By row and column.
    by.row <- sample(nrow(factors), nrow(factors)/2)
    by.col <- sample(ncol(factors), ncol(factors)/2)
    lem.alt <- lem[by.row, by.col]

    expect_identical(sampleFactors(lem.alt), factors[by.row,by.col])
    expect_identical(featureLoadings(lem.alt), loadings[,by.col])
    expect_identical(factorData(lem.alt), fdata[by.col,,drop=FALSE])

    # By row, with and without drop.
    keeper <- lem[1,]
    expect_identical(keeper, factors[1,])

    nodrop <- lem[1,,drop=FALSE]
    expect_identical(sampleFactors(nodrop), factors[1,,drop=FALSE])
    expect_identical(featureLoadings(nodrop), loadings)
    expect_identical(factorData(nodrop), fdata)

    # By column, with and without drop.
    keeper <- lem[,1]
    expect_identical(keeper, factors[,1])

    nodrop <- lem[,1,drop=FALSE]
    expect_identical(sampleFactors(nodrop), factors[,1,drop=FALSE])
    expect_identical(featureLoadings(nodrop), loadings[,1,drop=FALSE])
    expect_identical(factorData(nodrop), fdata[1,,drop=FALSE])

    # Throws errors correctly.
    expect_error(lem[nrow(lem)+1,], "subscript out of bounds", fixed=TRUE)
})

test_that("subsetting assignment works correctly", {
     # By row.
    src <- sample(nrow(lem), 10)
    dest <- sample(nrow(lem), 10)

    lem.alt <- lem
    lem.alt[dest,] <- lem[src,]

    ref <- factors
    ref[dest,] <- ref[src,]
    expect_identical(sampleFactors(lem.alt), ref)
    expect_identical(featureLoadings(lem.alt), loadings)
    expect_identical(factorData(lem.alt), fdata)

    # By column.
    src <- sample(ncol(lem), 5)
    dest <- sample(ncol(lem), 5)

    lem.alt <- lem
    lem.alt[,dest] <- lem[,src]

    ref_sf <- factors
    ref_sf[,dest] <- ref_sf[,src]
    expect_identical(sampleFactors(lem.alt), ref_sf)

    ref_fl <- loadings
    ref_fl[,dest] <- ref_fl[,src]
    expect_identical(featureLoadings(lem.alt), ref_fl)

    ref_fd <- fdata
    ref_fd[dest,] <- fdata[src,]
    expect_identical(factorData(lem.alt), ref_fd)

    # By row and column.
    dest.row <- sample(nrow(factors), nrow(factors)/2)
    src.row <- sample(nrow(factors), nrow(factors)/2)
    dest.col <- sample(ncol(factors), ncol(factors)/2)
    src.col <- sample(ncol(factors), ncol(factors)/2)
    
    lem.alt <- lem
    lem.alt[dest.row,dest.col] <- lem[src.row,src.col]

    ref_sf <- factors
    ref_sf[dest.row,dest.col] <- ref_sf[src.row,src.col]
    expect_identical(sampleFactors(lem.alt), ref_sf)

    ref_fl <- loadings
    ref_fl[,dest.col] <- ref_fl[,src.col]
    expect_identical(featureLoadings(lem.alt), ref_fl)

    ref_fd <- fdata
    ref_fd[dest.col,] <- fdata[src.col,]
    expect_identical(factorData(lem.alt), ref_fd)
})
