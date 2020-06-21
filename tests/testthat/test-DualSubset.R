# This tests the behavior of the DualSubset class.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-DualSubset.R")

set.seed(1000)
ds <- SingleCellExperiment:::DualSubset(rhits)
.hits <- SingleCellExperiment:::.get_hits

test_that("DualSubset subsetting works correctly", {
    REF_SUB <- function(x, i) {
        p <- alt <- SingleCellExperiment:::.get_hits(x)
        mcols(alt)[,1] <- seq_along(alt)
        mat <- SingleCellExperiment:::.hits2mat(alt)

        mat <- mat[i,i,drop=FALSE]
        p2 <- SingleCellExperiment:::.mat2hits(mat)
        mcols(p2) <- mcols(p)[mcols(p2)$x,,drop=FALSE]
        sort(p2)
    }

    keep <- 1:50
    sub <- ds[keep]
    expect_identical(length(sub), length(keep))
    expect_identical(.hits(sub), REF_SUB(ds, keep))
    expect_identical(
        length(.hits(sub)),
        sum(queryHits(rhits) %in% keep & subjectHits(rhits) %in% keep)
    )

    keep <- sample(length(ds), 50)
    sub <- ds[keep]
    expect_identical(length(sub), length(keep))
    expect_identical(.hits(sub), REF_SUB(ds, keep))
    expect_identical(
        length(.hits(sub)),
        sum(queryHits(rhits) %in% keep & subjectHits(rhits) %in% keep)
    )

    keep <- length(ds):1
    sub <- ds[keep]
    expect_identical(length(ds), length(sub))
    expect_identical(.hits(sub), REF_SUB(ds, keep))
    expect_identical(
        .hits(sub),
        sort(SelfHits(
            length(ds) - queryHits(rhits) + 1L,
            length(ds) - subjectHits(rhits) + 1L,
            nnode=length(ds),
            value=mcols(rhits)$value
        ))
    )
})

test_that("DualSubset subset replacement works correctly", {
    REF_SUB_REP <- function(x, i, value) {
        p <- alt <- SingleCellExperiment:::.get_hits(x)
        mcols(alt)[,1] <- seq_along(alt)
        mat <- SingleCellExperiment:::.hits2mat(alt)

        pv <- altv <- SingleCellExperiment:::.get_hits(value)
        mcols(altv)[,1] <- seq_along(altv)
        matv <- -SingleCellExperiment:::.hits2mat(altv)

        mat[i,i] <- matv
        p2 <- SingleCellExperiment:::.mat2hits(mat)
        index <- mcols(p2)$x
        use.left <- index > 0

        store <- mcols(p)[ifelse(use.left, index, 1),,drop=FALSE]
        store[!use.left,] <- mcols(pv)[-index[!use.left],,drop=FALSE]
        mcols(p2) <- store
        sort(p2)
    }

    swap1 <- 1:100
    swap2 <- 101:200

    ds2 <- ds
    ds2[swap1] <- ds[swap2]
    expect_identical(length(ds), length(ds2))
    expect_identical(.hits(ds2), REF_SUB_REP(ds, swap1, ds[swap2]))

    mod <- rhits
    mod <- mod[!(queryHits(mod) %in% swap1 & subjectHits(mod) %in% swap1)]
    sub <- .hits(ds[swap2])

    expect_identical(
        .hits(ds2),
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
    expect_identical(.hits(ds), .hits(ds2))
})

test_that("DualSubset concatenation works correctly", {
    combined <- c(ds, ds)
    expect_identical(length(combined), length(ds)*2L)
    expect_identical(
        .hits(combined),
        sort(SelfHits(
            c(queryHits(rhits), queryHits(rhits) + nnode(rhits)),
            c(subjectHits(rhits), subjectHits(rhits) + nnode(rhits)),
            nnode=nnode(rhits)*2L,
            value=c(mcols(rhits)$value, mcols(rhits)$value)
        ))
    )
})

test_that("SelfHits conversion to/from matrix works correctly", {
    p <- .hits(ds)

    mat <- SingleCellExperiment:::.hits2mat(p)
    expect_identical(nrow(mat), nnode(p))
    expect_identical(ncol(mat), nnode(p))
    expect_equal(sum(mat), sum(mcols(p)$value))

    roundtrip <- sort(SingleCellExperiment:::.mat2hits(mat))
    expect_identical(queryHits(roundtrip), queryHits(p))
    expect_identical(subjectHits(roundtrip), subjectHits(p))
    expect_identical(mcols(roundtrip)$x, mcols(p)$value)

    # Conversion still uses the first metadata field.
    mcols(p)$more <- "A"
    mat2 <- SingleCellExperiment:::.hits2mat(p)
    expect_identical(mat, mat2)

    mcols(p)$value <- NULL
    expect_error(SingleCellExperiment:::.hits2mat(p), "not supported")

    mcols(p)$more <- NULL
    mat2 <- SingleCellExperiment:::.hits2mat(p)
    expect_type(as.vector(mat2[0,0]), "logical")
    expect_equal(sum(mat2), length(p)) 
})

