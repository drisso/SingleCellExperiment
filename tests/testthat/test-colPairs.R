# Checks for proper functioning of the colPair methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-colPairs.R")

sce <- empty

chits <- sort(chits)
chits2 <- chits
mcols(chits2)$value <- mcols(chits2)$value * 2 
chits0 <- SelfHits(queryHits(chits), subjectHits(chits), nnode=nnode(chits) * 2L)

test_that("colPair getters/setters are functioning with character 'type'", {
    colPair(sce, "thing") <- chits
    expect_identical(colPair(sce, "thing"), chits)
    expect_identical(colPairNames(sce), "thing")

    colPair(sce, "stuff") <- chits2
    expect_identical(colPair(sce, "thing"), chits)
    expect_identical(colPair(sce, "stuff"), chits2)
    expect_identical(colPairNames(sce), c("thing", "stuff"))

    # Clearing values.
    colPair(sce, "thing") <- NULL
    expect_identical(colPair(sce, "stuff"), chits2)
    expect_identical(colPairNames(sce), "stuff")

    # In the absence of of redDim, create an unnamed one (like colPairs does)
    colPair(sce) <- NULL
    colPair(sce) <- chits
    expect_identical(colPairNames(sce), "unnamed1")

    # Checking for different errors.
    expect_error(colPair(sce, "thing"), "invalid subscript") 
    expect_error(colPair(sce, "DM") <- chits0, "number of nodes")
})

test_that("colPair getters/setters work with integer 'type'", {
    expect_error(colPair(sce), "no available entries") 
    expect_error(colPair(sce, 2), "invalid subscript 'type'") 

    expect_error(colPair(sce, 1) <- chits, "out of bounds")
    expect_error(colPair(sce, 2) <- chits, "out of bounds")

    # This gets a bit confusing as the order changes when earlier elements are wiped out.
    colPairs(sce) <- list(chits, chits2)
    expect_identical(colPair(sce), chits)
    expect_identical(colPair(sce, 2), chits2)
    expect_identical(colPairNames(sce), c("unnamed1", "unnamed2"))

    mult <- chits
    mcols(mult)$value <- mcols(mult)$value * 5
    colPair(sce, "thing") <- mult # chits is the second element.
    expect_identical(colPair(sce, 1), chits)
    expect_identical(colPair(sce, 2), chits2)
    expect_identical(colPair(sce, 3), mult)
    expect_identical(colPairNames(sce), c("unnamed1", "unnamed2", "thing"))

    colPair(sce, 1) <- NULL # chits2 becomes the first element now.
    expect_identical(colPair(sce), chits2)
    expect_identical(colPair(sce, 1), chits2)
    expect_identical(colPair(sce, 2), colPair(sce, "thing"))
    expect_identical(colPairNames(sce), c("unnamed2", "thing"))

    colPair(sce) <- NULL # 'mult' becomes the first element.
    expect_identical(colPair(sce), mult)
    expect_identical(colPairNames(sce), "thing")
    colPair(sce) <- chits2 # chits2 now overwrites the first element.
    expect_identical(colPair(sce, 1), chits2)
    expect_identical(colPairNames(sce), "thing")

    expect_error(colPair(sce, 5) <- chits, "out of bounds")
})

test_that("colPairs getters/setters are functioning", {
    colPairs(sce) <- list(thing=chits, TSNE=chits2)
    expect_identical(colPairs(sce), SimpleList(thing=chits, TSNE=chits2))

    expect_identical(colPairNames(sce), c("thing", "TSNE"))
    expect_identical(colPair(sce, "thing"), chits)
    expect_identical(colPair(sce, 1), chits)
    expect_identical(colPair(sce, "TSNE"), chits2)
    expect_identical(colPair(sce, 2), chits2)

    # Clearing via empty List.
    alt <- sce
    colPairs(alt) <- SimpleList()
    expect_identical(colPairs(alt), setNames(SimpleList(), character(0)))

    # Clearing via NULL.
    colPairs(sce) <- List(DM=chits)
    expect_identical(SimpleList(DM=chits), colPairs(sce))
    expect_identical(chits, colPair(sce))

    alt <- sce
    colPairs(alt) <- NULL
    expect_identical(colPairs(alt), setNames(SimpleList(), character(0)))

    # Setting with an unnamed list works.
    colPairs(sce) <- list(chits, chits2)
    expect_identical(colPairNames(sce), c("unnamed1", "unnamed2"))

    # Checking for errors.
    expect_error(colPairs(sce) <- list(chits, chits0), "number of nodes")
    expect_error(colPairs(sce) <- list(chits0, chits0), "number of nodes")
})

test_that("colPair getters/setters work with matrices", {
    mat <- SingleCellExperiment:::.hits2mat(chits)
    colPair(sce, "thing") <- mat
    expect_identical(colPair(sce, "thing", asSparse=TRUE), mat)

    colPairs(sce) <- list(thing=chits, TSNE=mat)
    everything <- colPairs(sce, asSparse=TRUE)
    expect_identical(everything[[1]], mat)

    colnames(mcols(chits)) <- "x"
    expect_identical(colPair(sce, 2), chits)
})

test_that("colPairNames getters/setters work correctly", {
    colPairs(sce) <- list(chits, chits2)
    expect_identical(colPairNames(sce), c("unnamed1", "unnamed2"))
    colPairs(sce) <- list(thing=chits, TSNE=chits2)
    expect_identical(colPairNames(sce), c("thing", "TSNE"))

    # Directly setting.
    colPairNames(sce) <- c("A", "B")
    expect_identical(colPairNames(sce), c("A", "B"))

    # When wiped.
    colPairs(sce) <- NULL
    expect_identical(colPairNames(sce), character(0))

    expect_error(colPairNames(sce) <- c("A", "B"), "more column names")
})
