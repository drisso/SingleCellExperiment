# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("test-lem-combine.R")
context("LEM combine")

set.seed(1000)
ncells <- 100

factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
fdata <- DataFrame(WHEE=sample(LETTERS, ncol(factors)))
lem <- LinearEmbeddingMatrix(factors, loadings, fdata)

test_that("rbind works correctly", {
    shuffled <- sample(nrow(lem))
    lem.alt <- lem[shuffled,]
    samish <- rbind(lem.alt)
    expect_identical(sampleFactors(lem.alt), sampleFactors(samish))
    expect_identical(featureLoadings(lem.alt), featureLoadings(samish))
    expect_identical(factorData(lem.alt), factorData(samish))

    lem2 <- rbind(lem, lem.alt)
    expect_identical(sampleFactors(lem2, withDimnames=FALSE),
                     rbind(sampleFactors(lem, withDimnames=FALSE),
                           sampleFactors(lem.alt, withDimnames=FALSE)))
    expect_identical(featureLoadings(lem2), featureLoadings(lem))
    expect_identical(factorData(lem2), factorData(lem))
    expect_identical(rownames(lem2), NULL)
    expect_identical(colnames(lem2), NULL)

    # Works correctly with names.
    unnamed <- lem
    rownames(lem) <- paste0("CELL", seq_len(nrow(lem)))
    colnames(lem) <- paste0("FACTOR", seq_len(ncol(lem)))

    lem3 <- rbind(lem, lem[shuffled,])
    expect_identical(rownames(lem3), c(rownames(lem), rownames(lem)[shuffled]))
    expect_identical(colnames(lem3), colnames(lem))
    expect_equivalent(sampleFactors(lem3), sampleFactors(lem2))

    lem4 <- rbind(lem, unnamed)
    expect_identical(rownames(lem4), NULL)
    expect_identical(colnames(lem4), colnames(lem))
    lem4 <- rbind(unnamed, lem)
    expect_identical(rownames(lem4), NULL)
    expect_identical(colnames(lem4), NULL)

    # Throws errors correctly.
    lem.alt <- lem[,1:2]
    expect_error(rbind(lem, lem.alt), "number of columns of matrices must match")
})

test_that("cbind works correctly", {
    shuffled <- sample(ncol(factors))
    lem.alt <- lem[,shuffled]
    samish <- cbind(lem.alt)
    expect_identical(sampleFactors(lem.alt), sampleFactors(samish))
    expect_identical(featureLoadings(lem.alt), featureLoadings(samish))
    expect_identical(factorData(lem.alt), factorData(samish))

    lem2 <- cbind(lem, lem.alt)
    expect_identical(sampleFactors(lem2, withDimnames=FALSE),
                     cbind(sampleFactors(lem, withDimnames=FALSE),
                           sampleFactors(lem.alt, withDimnames=FALSE)))
    expect_identical(featureLoadings(lem2), cbind(featureLoadings(lem), featureLoadings(lem.alt)))
    expect_identical(factorData(lem2), rbind(factorData(lem), factorData(lem.alt)))
    expect_identical(rownames(lem2), NULL)
    expect_identical(colnames(lem2), NULL)

    # Works correctly with names.
    unnamed <- lem
    rownames(lem) <- paste0("CELL", seq_len(nrow(lem)))
    colnames(lem) <- paste0("FACTOR", seq_len(ncol(lem)))

    lem3<- cbind(lem, lem[,shuffled])
    expect_identical(rownames(lem3), rownames(lem))
    expect_identical(colnames(lem3), c(colnames(lem), colnames(lem)[shuffled]))
    expect_equivalent(sampleFactors(lem3), sampleFactors(lem2))

    lem4 <- cbind(unnamed, lem)
    expect_identical(rownames(lem4), NULL)
    expect_identical(colnames(lem4), NULL)
    lem4 <- cbind(lem, unnamed)
    expect_identical(rownames(lem4), rownames(lem))
    expect_identical(colnames(lem4), NULL)

    # Throws errors correctly.
    lem.alt <- lem[1:5,]
    expect_error(cbind(lem, lem.alt), "number of rows")

    lem.alt <- lem
    factorData(lem.alt)$WHEE <- NULL
    expect_error(cbind(lem, lem.alt), "column names")
})
