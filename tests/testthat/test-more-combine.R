# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); source("setup.R"); source("test-more-combine.R")

sce <- loaded
sce2 <- loaded[,1:50]

test_that("basic combineCols works correctly", {
    copy <- sce2
    logcounts(copy) <- NULL

    out <- combineCols(sce, copy, use.names=FALSE)
    expect_identical(ncol(out), ncol(sce) + ncol(copy))
    expect_identical(as.matrix(logcounts(out)[,1:ncol(sce)]), as.matrix(logcounts(sce)))
    expect_true(all(is.na(as.matrix(logcounts(out)[,ncol(sce) + 1:ncol(copy)]))))

    rownames(sce2) <- rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
    copy <- sce2[10:20,]
    out <- combineCols(sce, copy, use.names=TRUE)
    expect_identical(ncol(out), ncol(sce) + ncol(copy))
    expect_identical(as.matrix(counts(out)[,1:ncol(sce)]), as.matrix(counts(sce)))
    expect_identical(as.matrix(counts(out)[10:20,ncol(sce) + 1:ncol(copy)]), as.matrix(counts(copy)))
    expect_true(all(is.na(counts(out)[-(10:20),ncol(sce) + 1:ncol(copy)])))
})

test_that("combineCols works correctly with reducedDims", {
    # Handles absent entries.
    copy <- sce2
    reducedDim(copy, "PCA") <- NULL
    reducedDim(copy, "TSNE") <- reducedDim(copy, "TSNE") * 2
    out <- combineCols(sce, copy, use.names=FALSE)

    expect_identical(reducedDim(out, "TSNE"), rbind(reducedDim(sce, "TSNE"), reducedDim(copy, "TSNE")))
    expect_identical(reducedDim(out, "PCA")[1:ncol(sce),], reducedDim(sce, "PCA"))
    expect_true(all(is.na(reducedDim(out, "PCA")[ncol(sce) + 1:ncol(copy),])))

    # Robust to incompatible dimensions.
    copy <- sce2
    reducedDim(copy, "PCA") <- reducedDim(copy, "PCA")[,1:2]
    expect_warning(out <- combineCols(sce, copy, use.names=FALSE), "PCA")
    expect_false("PCA" %in% reducedDimNames(out))
    expect_identical(reducedDim(out, "TSNE"), rbind(reducedDim(sce, "TSNE"), reducedDim(copy, "TSNE")))
})

test_that("combineCols works correctly with altExps", {
    # Handles absent entries.
    copy1 <- sce
    altExp(copy1, "Protein") <- NULL
    copy2 <- sce2
    assay(altExp(copy2, "Spike"), "counts") <- NULL
    out <- combineCols(copy1, copy2, use.names=FALSE)

    prot <- altExp(out, "Protein")
    expect_true(all(is.na(as.matrix(assay(prot, "logcounts")[,1:ncol(copy1)]))))
    expect_identical(as.matrix(assay(prot, "logcounts")[,ncol(copy1) + 1:ncol(copy2)]), as.matrix(assay(altExp(copy2, "Protein"), "logcounts")))

    spike <- altExp(out, "Spike")
    expect_identical(as.matrix(assay(spike, "counts")[,1:ncol(copy1)]), as.matrix(assay(altExp(copy1, "Spike"), "counts")))
    expect_true(all(is.na(as.matrix(assay(spike, "counts")[,ncol(copy1) + 1:ncol(copy2)]))))

    # Robust to incompatible dimensions.
    copy1 <- sce
    altExp(copy1, "Protein") <- altExp(copy1, "Protein")[1:2,]
    expect_warning(out <- combineCols(copy1, copy2, use.names=FALSE), "Protein")
    expect_false("Protein" %in% altExpNames(out))

    # Use.names=TRUE propagates.
    rownames(copy1) <- rownames(copy2) <- sprintf("GENE_%i", seq_len(nrow(copy1)))
    out <- combineCols(copy1, copy2)
    prot <- altExp(out, "Protein")
    vals <- assay(prot, "counts")
    expect_identical(as.matrix(vals[1:2,1:ncol(copy1)]), as.matrix(assay(altExp(copy1, "Protein"), "counts")))
    expect_true(all(is.na(vals[-(1:2),1:ncol(copy1)])))
})

test_that("combineCols works correctly with colPairs", {
    # Handles absent entries.
    copy1 <- sce
    colPairs(copy1) <- NULL

    out <- combineCols(copy1, sce2, use.names=FALSE)
    expect_identical(as.matrix(colPair(out)) - ncol(copy1), as.matrix(colPair(sce2)))

    out <- combineCols(sce, copy1, use.names=FALSE)
    expect_identical(as.data.frame(colPair(out)), as.data.frame(colPair(sce)))
})

test_that("combineCols works correctly with rowPairs", {
    # Handles absent entries.
    copy1 <- sce
    rowPairs(copy1) <- NULL
    out <- combineCols(copy1, sce2, use.names=FALSE)
    expect_identical(rowPairs(out), rowPairs(sce2))

    # Handles mismatching entries.
    copy1 <- sce
    mcols(rowPair(copy1))$value <- runif(100)
    expect_warning(out <- combineCols(copy1, sce2, use.names=FALSE), "unnamed1")
    expect_identical(rowPairs(out), rowPairs(copy1))

    # Merges non-contradictory entries.
    copy1 <- sce
    copy2 <- sce2
    rownames(copy1) <- rownames(copy2) <- sprintf("GENE_%i", seq_len(nrow(copy1)))
    copy1 <- copy1[1:(nrow(sce)/2),]
    copy2 <- copy2[nrow(sce)/2+1:(nrow(sce)/2),]

    expect_warning(out <- combineCols(copy1, copy2), NA)
    expect_true(length(rowPair(out)) > 0)
    expect_true(length(rowPair(out)) < length(rowPair(sce))) # as links spanning the split are lost.
})
