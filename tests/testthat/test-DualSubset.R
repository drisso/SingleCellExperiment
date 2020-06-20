# This tests the behavior of the DualSubset class.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-DualSubset.R")

set.seed(1000)
ds <- SingleCellExperiment:::DualSubset(rhits)

test_that("DualSubset subsetting works correctly", {
    keep <- 1:50
    expect_identical(length(ds[keep]), length(keep))
    expect_identical(
        length(SingleCellExperiment:::.get_hits(ds[keep])),
        sum(queryHits(rhits) %in% keep & subjectHits(rhits) %in% keep)
    )

    keep <- sample(length(ds), 50)
    expect_identical(length(ds[keep]), length(keep))
    expect_identical(
        length(SingleCellExperiment:::.get_hits(ds[keep])),
        sum(queryHits(rhits) %in% keep & subjectHits(rhits) %in% keep)
    )

    keep <- length(ds):1
    expect_identical(length(ds), length(ds[keep]))
    expect_identical(
        sort(SingleCellExperiment:::.get_hits(ds[keep])),
        sort(SelfHits(
            length(ds) - queryHits(rhits) + 1L,
            length(ds) - subjectHits(rhits) + 1L,
            nnode=length(ds),
            value=mcols(rhits)$value
        ))
    )
})

test_that("DualSubset subset replacement works correctly", {
    swap1 <- 1:100
    swap2 <- 101:200

    ds2 <- ds
    ds2[swap1] <- ds[swap2]
    expect_identical(length(ds), length(ds2))

    mod <- rhits
    mod <- mod[!(queryHits(mod) %in% swap1 & subjectHits(mod) %in% swap1)]
    sub <- SingleCellExperiment:::.get_hits(ds[swap2])

    expect_identical(
        sort(SingleCellExperiment:::.get_hits(ds2)),
        sort(SelfHits(
            c(queryHits(mod), queryHits(sub)),
            c(subjectHits(mod), subjectHits(sub)),
            nnode=length(ds),
            value=c(mcols(mod)$value, mcols(sub)$value)
        ))
    )

    # This should be a no-op.
    ds2 <- ds
    ds2[swap1] <- ds2[swap1]
    expect_identical(ds, ds2)
})

test_that("DualSubset concatenation works correctly", {
    combined <- c(ds, ds)
    expect_identical(length(combined), length(ds)*2L)
    expect_identical(
        sort(SingleCellExperiment:::.get_hits(combined)),
        sort(SelfHits(
            c(queryHits(rhits), queryHits(rhits) + nnode(rhits)),
            c(subjectHits(rhits), subjectHits(rhits) + nnode(rhits)),
            nnode=nnode(rhits)*2L,
            value=c(mcols(rhits)$value, mcols(rhits)$value)
        ))
    )
})
