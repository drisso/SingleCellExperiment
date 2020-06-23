# Checks for proper functioning of the rowPair methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-rowPairs.R")

sce <- empty

rhits <- sort(rhits)
rhits2 <- rhits
mcols(rhits2)$value <- mcols(rhits2)$value * 2 
rhits0 <- SelfHits(queryHits(rhits), subjectHits(rhits), nnode=nnode(rhits) * 2L)

test_that("rowPair getters/setters are functioning with character 'type'", {
    rowPair(sce, "thing") <- rhits
    expect_identical(rowPair(sce, "thing"), rhits)
    expect_identical(rowPairNames(sce), "thing")

    rowPair(sce, "stuff") <- rhits2
    expect_identical(rowPair(sce, "thing"), rhits)
    expect_identical(rowPair(sce, "stuff"), rhits2)
    expect_identical(rowPairNames(sce), c("thing", "stuff"))

    # Clearing values.
    rowPair(sce, "thing") <- NULL
    expect_identical(rowPair(sce, "stuff"), rhits2)
    expect_identical(rowPairNames(sce), "stuff")

    # In the absence of of redDim, create an unnamed one (like rowPairs does)
    rowPair(sce) <- NULL
    rowPair(sce) <- rhits
    expect_identical(rowPairNames(sce), "unnamed1")

    # Checking for different errors.
    expect_error(rowPair(sce, "thing"), "invalid subscript") 
    expect_error(rowPair(sce, "DM") <- rhits0, "number of nodes")
})

test_that("rowPair getters/setters work with integer 'type'", {
    expect_error(rowPair(sce), "no available entries") 
    expect_error(rowPair(sce, 2), "invalid subscript 'type'") 

    expect_error(rowPair(sce, 1) <- rhits, "out of bounds")
    expect_error(rowPair(sce, 2) <- rhits, "out of bounds")

    # This gets a bit confusing as the order changes when earlier elements are wiped out.
    rowPairs(sce) <- list(rhits, rhits2)
    expect_identical(rowPair(sce), rhits)
    expect_identical(rowPair(sce, 2), rhits2)
    expect_identical(rowPairNames(sce), c("unnamed1", "unnamed2"))

    mult <- rhits
    mcols(mult)$value <- mcols(mult)$value * 5
    rowPair(sce, "thing") <- mult # rhits is the second element.
    expect_identical(rowPair(sce, 1), rhits)
    expect_identical(rowPair(sce, 2), rhits2)
    expect_identical(rowPair(sce, 3), mult)
    expect_identical(rowPairNames(sce), c("unnamed1", "unnamed2", "thing"))

    rowPair(sce, 1) <- NULL # rhits2 becomes the first element now.
    expect_identical(rowPair(sce), rhits2)
    expect_identical(rowPair(sce, 1), rhits2)
    expect_identical(rowPair(sce, 2), rowPair(sce, "thing"))
    expect_identical(rowPairNames(sce), c("unnamed2", "thing"))

    rowPair(sce) <- NULL # 'mult' becomes the first element.
    expect_identical(rowPair(sce), mult)
    expect_identical(rowPairNames(sce), "thing")
    rowPair(sce) <- rhits2 # rhits2 now overwrites the first element.
    expect_identical(rowPair(sce, 1), rhits2)
    expect_identical(rowPairNames(sce), "thing")

    expect_error(rowPair(sce, 5) <- rhits, "out of bounds")
})

test_that("rowPairs getters/setters are functioning", {
    rowPairs(sce) <- list(thing=rhits, TSNE=rhits2)
    expect_identical(rowPairNames(sce), c("thing", "TSNE"))
    expect_identical(rowPair(sce, "thing"), rhits)
    expect_identical(rowPair(sce, 1), rhits)
    expect_identical(rowPair(sce, "TSNE"), rhits2)
    expect_identical(rowPair(sce, 2), rhits2)

    # Clearing via empty List.
    alt <- sce
    rowPairs(alt) <- SimpleList()
    expect_identical(rowPairs(alt), setNames(SimpleList(), character(0)))

    # Clearing via NULL.
    rowPairs(sce) <- List(DM=rhits)
    expect_identical(SimpleList(DM=rhits), rowPairs(sce))
    expect_identical(rhits, rowPair(sce))

    alt <- sce
    rowPairs(alt) <- NULL
    expect_identical(rowPairs(alt), setNames(SimpleList(), character(0)))

    # Setting with an unnamed list works.
    rowPairs(sce) <- list(rhits, rhits2)
    expect_identical(rowPairNames(sce), c("unnamed1", "unnamed2"))

    # Checking for errors.
    expect_error(rowPairs(sce) <- list(rhits, rhits0), "number of nodes")
    expect_error(rowPairs(sce) <- list(rhits0, rhits0), "number of nodes")
})

test_that("rowPair getters/setters work with matrices", {
    mat <- SingleCellExperiment:::.hits2mat(rhits)
    rowPair(sce, "thing") <- mat
    expect_identical(rowPair(sce, "thing", asSparse=TRUE), mat)

    rowPairs(sce) <- list(thing=rhits, TSNE=mat)
    everything <- rowPairs(sce, asSparse=TRUE)
    expect_identical(everything[[1]], mat)

    colnames(mcols(rhits)) <- "x"
    expect_identical(rowPair(sce, 2), rhits)
})

test_that("rowPairNames getters/setters work correctly", {
    rowPairs(sce) <- list(rhits, rhits2)
    expect_identical(rowPairNames(sce), c("unnamed1", "unnamed2"))
    rowPairs(sce) <- list(thing=rhits, TSNE=rhits2)
    expect_identical(rowPairNames(sce), c("thing", "TSNE"))

    # Directly setting.
    rowPairNames(sce) <- c("A", "B")
    expect_identical(rowPairNames(sce), c("A", "B"))

    # When wiped.
    rowPairs(sce) <- NULL
    expect_identical(rowPairNames(sce), character(0))

    expect_error(rowPairNames(sce) <- c("A", "B"), "more column names")
})
