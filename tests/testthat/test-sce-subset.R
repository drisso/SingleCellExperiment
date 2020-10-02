# Checks the subsetting methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-subset.R")

sce <- loaded
rownames(sce) <- paste0("Gene", seq_len(nrow(v)))
colnames(sce) <- paste0("Cell", seq_len(ncells))

test_that("subsetting by row works correctly", {
    int_elementMetadata(sce)$indicator <- seq_len(nrow(sce))

    for (i in 1:3) {
        if (i==1L) {
            by.row <- sample(nrow(v), 20)
            sub.sce <- sce[by.row,]
            expect_identical(int_elementMetadata(sub.sce)$indicator, by.row)
        } else if (i==2L) {
            by.row <- rbinom(nrow(v), 1, 0.2)==1
            sub.sce <- sce[by.row,]
            expect_identical(int_elementMetadata(sub.sce)$indicator, which(by.row))
        } else if (i==3L) {
            by.row <- rownames(sce)[sample(nrow(v), 100)]
            sub.sce <- sce[by.row,]
            expect_identical(int_elementMetadata(sub.sce)$indicator, match(by.row, rownames(sce)))
        }
        ind <- int_elementMetadata(sub.sce)$indicator

        expect_identical(assay(sce)[ind,,drop=FALSE], assay(sub.sce)) # check SE elements are subsetted.
        expect_identical(rowData(sce)[ind,], rowData(sub.sce))
        expect_identical(rowPair(sub.sce), 
            SingleCellExperiment:::.get_hits(SingleCellExperiment:::DualSubset(rowPair(sce))[ind]))

        # Unchanged elements:
        expect_identical(sizeFactors(sub.sce), sizeFactors(sce))
        expect_identical(reducedDims(sub.sce), reducedDims(sce))
        expect_identical(int_colData(sub.sce), int_colData(sce))
        expect_identical(objectVersion(sub.sce), objectVersion(sce))
    }

    expect_error(sce[nrow(sce)+1,], "subscript contains out-of-bounds indices", fixed=TRUE)
    expect_error(sce["A",], "index out of bounds: A")
})

test_that("subsetting by column works correctly", {
    colData(sce)$indicator <- seq_len(ncells)

    for (j in 1:3) {
        if (j==1L) {
            by.col <- sample(ncells, 20)
            sub.sce <- sce[,by.col]
            expect_identical(colData(sub.sce)$indicator, by.col)
        } else if (j==2L) {
            by.col <- rbinom(ncells, 1, 0.2)==1
            sub.sce <- sce[,by.col]
            expect_identical(colData(sub.sce)$indicator, which(by.col))
        } else if (j==3L) {
            by.col <- colnames(sce)[sample(ncells, 50)]
            sub.sce <- sce[,by.col]
            expect_identical(colData(sub.sce)$indicator, match(by.col, colnames(sce)))
        }
        ind <- colData(sub.sce)$indicator

        expect_identical(assay(sce)[,ind,drop=FALSE], assay(sub.sce)) # check SE elements are subsetted.
        expect_identical(sizeFactors(sub.sce), sizeFactors(sce)[ind])

        expect_identical(reducedDim(sub.sce, "PCA", withDimnames=FALSE), d1[ind,,drop=FALSE])
        expect_identical(reducedDim(sub.sce, "TSNE", withDimnames=FALSE), d2[ind,,drop=FALSE])

        expect_identical(colPair(sub.sce),
            SingleCellExperiment:::.get_hits(SingleCellExperiment:::DualSubset(colPair(sce))[ind]))

        # Unchanged elements:
        expect_identical(int_elementMetadata(sub.sce), int_elementMetadata(sce))
        expect_identical(rowData(sub.sce), rowData(sce))
        expect_identical(objectVersion(sub.sce), objectVersion(sce))
    }

    expect_error(sce[,ncells+1], "subscript contains out-of-bounds indices", fixed=TRUE)
    expect_error(sce[,"A"], "index out of bounds: A")
})

test_that("subset replacement by row works correctly for basic cases", {
    sce.alt <- sce
    rownames(sce.alt) <- paste0(rownames(sce), "x")
    int_metadata(sce.alt)$whee <- 1

    scex <- sce.alt
    to <- 1:10
    from <- 21:30
    scex[to,] <- sce[from,]

    expect_identical(assay(scex)[to,,drop=FALSE], assay(sce)[from,,drop=FALSE])
    expect_equivalent(assay(scex)[-to,,drop=FALSE], assay(sce)[-to,,drop=FALSE])

    sub <- SingleCellExperiment:::DualSubset(rowPair(scex)) 
    full <- SingleCellExperiment:::DualSubset(rowPair(sce))
    expect_identical(sub[to], full[from])
    expect_identical(sub[-to], full[-to])

    # Unchanged elements, to name a few. 
    expect_identical(int_metadata(scex), int_metadata(sce))
    expect_identical(int_colData(scex), int_colData(sce))
    expect_identical(colData(scex), colData(sce))

    # Again for character
    scex2 <- sce.alt 
    to <- rownames(scex2)[1:10]
    from <- rownames(sce)[21:30]
    scex2[to,] <- sce[from,]
    expect_equal(scex, scex2)
})

test_that("subset replacement by row handles internal fields correctly", {
    to <- 1:10
    from <- 21:30

    # Handles mismatch.
    scex2 <- sce
    int_elementMetadata(scex2)$ERCC <- seq_len(nrow(scex2))
    expect_error(scex2[to,] <- sce[from,], "'int_elementMetadata'")
})

test_that("subset replacement by column works correctly for basic cases", {
    sce.alt <- sce
    colnames(sce.alt) <- paste0(colnames(sce), "x")

    scex <- sce.alt
    to <- 1:10
    from <- 21:30
    scex[,to] <- sce[,from]

    expect_identical(assay(scex)[,to,drop=FALSE], assay(sce)[,from,drop=FALSE])
    expect_equivalent(assay(scex)[,-to,drop=FALSE], assay(sce)[,-to,drop=FALSE])
    expect_identical(sizeFactors(scex)[to], sizeFactors(sce)[from])
    expect_equivalent(sizeFactors(scex)[-to], sizeFactors(sce)[-to])

    # Unchanged elements.
    expect_identical(int_elementMetadata(scex), int_elementMetadata(sce))
    expect_identical(rowData(scex), rowData(sce))
    expect_identical(int_metadata(scex), int_metadata(sce))
})

test_that("subset replacement by column handles internal fields", {
    to <- 1:10
    from <- 21:30

    # Throws an error upon mismatch.
    scex2 <- sce
    reducedDims(scex2) <- NULL
    expect_error(scex2[,to] <- sce[,from], "'int_colData'")
})

test_that("subset replacement by both rows and columns work correctly", {
    # Wholesale replacement!
    sce.alt <- sce
    rownames(sce.alt) <- paste0(rownames(sce), "x")
    colnames(sce.alt) <- paste0(colnames(sce), "x")
    int_metadata(sce.alt)$whee <- 1

    scex <- sce.alt
    scex[] <- sce
    expect_equal(scex, sce)

    # Partial replacement.
    to <- 1:10
    from <- 21:30

    scex <- sce.alt
    scex[to,to] <- sce[from,from]
    
    ref <- sce.alt
    ref[to,] <- sce[from,]
    expect_identical(int_elementMetadata(ref), int_elementMetadata(scex))

    ref <- sce.alt
    ref[,to] <- sce[,from]
    expect_identical(int_colData(ref), int_colData(scex))
})

test_that("S4Vectors subsetting works correctly", {
    out <- extractROWS(sce, 1:10)
    expect_identical(out, sce[1:10,])

    set.seed(100)
    f <- sample(10, nrow(sce), replace=TRUE)
    out <- split(sce, f)
    expect_identical(out[["1"]], sce[f==1,])
    expect_identical(out[["5"]], sce[f==5,])
    expect_identical(out[["8"]], sce[f==8,])
})
