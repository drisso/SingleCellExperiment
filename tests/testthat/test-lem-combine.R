# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("test-lem-combine.R")

set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
fdata <- DataFrame(WHEE=sample(LETTERS, ncol(factors)))
lem <- LinearEmbeddingMatrix(factors, loadings, fdata)

test_that("rbind works correctly", {
    shuffled <- sample(nrow(lem))
    lem.alt <- lem[shuffled,]
    expect_identical(lem.alt, rbind(lem.alt))

    lem2 <- rbind(lem, lem.alt)
    expect_identical(sampleFactors(lem2), rbind(sampleFactors(lem), sampleFactors(lem.alt)))
    expect_identical(featureLoadings(lem2), featureLoadings(lem))
    expect_identical(factorData(lem2), factorData(lem))

    # Throws errors correctly.
    lem.alt <- lem[,1:2]
    expect_error(rbind(lem, lem.alt), "number of columns of matrices must match")
})

test_that("cbind works correctly", {
    shuffled <- sample(ncol(factors))
    lem.alt <- lem[,shuffled]
    expect_identical(lem.alt, cbind(lem.alt))
    
    lem2 <- cbind(lem, lem.alt)
    expect_identical(sampleFactors(lem2), cbind(sampleFactors(lem), sampleFactors(lem.alt)))
    expect_identical(featureLoadings(lem2), cbind(featureLoadings(lem), featureLoadings(lem.alt)))
    expect_identical(factorData(lem2), rbind(factorData(lem), factorData(lem.alt)))

    # Throws errors correctly.
    lem.alt <- lem[1:5,]
    expect_error(cbind(lem, lem.alt), "number of rows")

    lem.alt <- lem
    factorData(lem.alt)$WHEE <- NULL
    expect_error(cbind(lem, lem.alt), "number of columns")
})
