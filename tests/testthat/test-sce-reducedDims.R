# Checks for proper functioning of the reducedDim methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-sce-reducedDims.R")

sce <- empty

test_that("reducedDim getters/setters are functioning with character 'type'", {
    reducedDim(sce, "PCA") <- d1
    expect_identical(reducedDim(sce, "PCA"), d1)
    expect_identical(reducedDims(sce), SimpleList(PCA=d1))
    expect_identical(reducedDimNames(sce), "PCA")

    reducedDim(sce, "tSNE") <- d2
    expect_identical(reducedDim(sce, "tSNE"), d2)
    expect_identical(reducedDims(sce), SimpleList(PCA=d1, tSNE=d2))
    expect_identical(reducedDimNames(sce), c("PCA", "tSNE"))

    # Clearing values.
    reducedDim(sce, "PCA") <- NULL
    expect_identical(reducedDim(sce, "tSNE"), d2)
    expect_identical(reducedDim(sce), d2)
    expect_identical(reducedDims(sce), SimpleList(tSNE=d2))
    expect_identical(reducedDimNames(sce), "tSNE")

    # Checking for different errors.
    expect_error(reducedDim(sce, "PCA"), "invalid names")
    expect_error(reducedDim(sce, 2), "out-of-bounds indices")
    expect_error(reducedDim(sce, "DM") <- d1[1:10,], "different number of rows")
    expect_error(reducedDim(sce, "DM") <- "huh", "different number of rows")
})

test_that("reducedDims getters/setters are functioning", {
    reducedDims(sce) <- list(PCA=d1, TSNE=d2)
    expect_identical(reducedDimNames(sce), c("PCA", "TSNE"))
    expect_identical(reducedDim(sce, "PCA"), d1)
    expect_identical(reducedDim(sce, 1), d1)
    expect_identical(reducedDim(sce, "TSNE"), d2)
    expect_identical(reducedDim(sce, 2), d2)

    # Clearing via empty List.
    alt <- sce
    reducedDims(alt) <- SimpleList()
    expect_identical(reducedDims(alt), setNames(SimpleList(), character(0)))

    # Clearing via NULL.
    reducedDims(sce) <- SimpleList(DM=d1)
    expect_identical(SimpleList(DM=d1), reducedDims(sce))
    expect_identical(d1, reducedDim(sce))

    alt <- sce
    reducedDims(alt) <- NULL
    expect_identical(reducedDims(alt), setNames(SimpleList(), character(0)))

    # Setting with an unnamed list works. 
    reducedDims(sce) <- list(d1, d2)
    expect_identical(length(reducedDimNames(sce)), 2L)

    # Checking for errors.
    expect_error(reducedDims(sce) <- list(d1, d2[1:10,]), "do not have the correct number of rows")
    expect_error(reducedDims(sce) <- list(d1[1:10,], d2[1:10,]), "do not have the correct number of rows")
})

test_that("getters/setters respond to dimnames", {
    named <- sce
    colnames(named) <- paste0("Cell", sample(ncol(named)))

    reducedDim(named, "PCA") <- d1
    reducedDim(named, "tSNE") <- d2
    expect_identical(rownames(reducedDim(named)), colnames(named))
    expect_identical(rownames(reducedDim(named, 2)), colnames(named))
    expect_identical(rownames(reducedDim(named, withDimnames=FALSE)), NULL)

    out <- reducedDims(named)
    expect_identical(rownames(out[[1]]), colnames(named))
    expect_identical(rownames(out[[2]]), colnames(named))
    out <- reducedDims(named, withDimnames=FALSE)
    expect_identical(rownames(out[[1]]), NULL)
    expect_identical(rownames(out[[2]]), NULL)
})

test_that("reducedDim getters/setters work with numeric indices", {
    # This gets a bit confusing as the order changes when earlier elements are wiped out.
    reducedDim(sce, 1) <- d1
    expect_identical(reducedDim(sce), d1)
    expect_identical(reducedDim(sce, 1), d1)
    expect_identical(reducedDimNames(sce), "")
    reducedDim(sce, 2) <- d2
    expect_identical(reducedDim(sce), d1)
    expect_identical(reducedDim(sce, 2), d2)
    expect_identical(reducedDimNames(sce), character(2))

    mult <- d1 * 5
    reducedDim(sce, "PCA") <- mult # d1 is the second element.
    expect_identical(reducedDim(sce, 1), d1)
    expect_identical(reducedDim(sce, 2), d2)
    expect_identical(reducedDim(sce, 3), mult)
    expect_identical(reducedDimNames(sce), c("", "", "PCA"))

    reducedDim(sce, 1) <- NULL # d2 becomes the first element now.
    expect_identical(reducedDim(sce), d2)
    expect_identical(reducedDim(sce, 1), d2)
    expect_identical(reducedDim(sce, 2), reducedDim(sce, "PCA"))
    expect_identical(reducedDimNames(sce), c("", "PCA"))

    reducedDim(sce) <- NULL # 'mult' becomes the first element.
    expect_identical(reducedDim(sce), mult)
    expect_identical(reducedDimNames(sce), "PCA")
    reducedDim(sce) <- d2 # d2 now overwrites the first element.
    expect_identical(reducedDim(sce, 1), d2)
    expect_identical(reducedDimNames(sce), "PCA")

    expect_error(reducedDim(sce, 5) <- d1, "subscript out of bounds")
})

test_that("reducedDimNames getters/setters work correctly", {
    reducedDims(sce) <- list(d1, d2)
    expect_identical(reducedDimNames(sce), character(2))
    reducedDims(sce) <- list(PCA=d1, TSNE=d2)
    expect_identical(reducedDimNames(sce), c("PCA", "TSNE"))

    # Directly setting.
    reducedDimNames(sce) <- c("A", "B")
    expect_identical(reducedDimNames(sce), c("A", "B"))

    # When wiped.
    reducedDims(sce) <- NULL
    expect_identical(reducedDimNames(sce), character(0))

    expect_error(reducedDimNames(sce) <- c("A", "B"), "more column names")
})
