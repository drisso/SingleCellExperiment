# Tests the unsplitAltExps function.
# library(testthat); library(SingleCellExperiment); source("test-unsplit-altexp.R")

#########################
#########################

# Assembling an appropriately stripped down dummy object.

set.seed(1000)
ncells <- 100
u <- matrix(rpois(20000, 5), ncol=ncells)
sce <- SingleCellExperiment(list(counts=u))
rownames(sce) <- sprintf("GENE_%i", seq_len(nrow(sce)))
sce$Cell_Cycle <- sample(LETTERS, ncol(sce), replace=TRUE)

se1 <- SingleCellExperiment(list(counts=matrix(rpois(2000, 5), ncol=ncells)))
rowData(se1)$stuff <- sample(LETTERS, nrow(se1), replace=TRUE)
rownames(se1) <- sprintf("SPIKE_%i", seq_len(nrow(se1)))
altExp(sce, "Spike") <- se1

se2 <- SummarizedExperiment(list(counts=matrix(rpois(800, 5), ncol=ncells)))
rowData(se2)$blah <- sample(letters, nrow(se2), replace=TRUE)
rownames(se2) <- sprintf("TAG_%i", seq_len(nrow(se2)))
se2$something <- sample(letters, ncol(sce), replace=TRUE)
altExp(sce, "Protein") <- se2

#########################
#########################

is.first <- seq_len(nrow(sce))
is.spike <- nrow(sce) + seq_len(nrow(altExp(sce)))
is.stuff <- nrow(sce) + nrow(altExp(sce)) + + seq_len(nrow(altExp(sce, 2)))

expected_names <- c(
    rownames(sce), 
    paste0("Spike.", rownames(altExp(sce))),
    paste0("Protein.", rownames(altExp(sce, 2)))
)

test_that("unsplitAltExps works as expected for assays", {
    output <- unsplitAltExps(sce, prefix.row=FALSE)
    expect_identical(as.matrix(assay(output)), rbind(assay(sce), assay(altExp(sce)), assay(altExp(sce, 2))))

    output <- unsplitAltExps(sce)
    expect_identical(expected_names, rownames(output))

    # A trickier situation with an assay only available in one of the matrices.
    alt <- sce
    logcounts(alt) <- log2(counts(sce) + 1)
    output <- unsplitAltExps(alt)
    expect_identical(as.matrix(logcounts(alt)), as.matrix(logcounts(output)[is.first,]))
    expect_true(all(is.na(logcounts(output)[is.spike,])))
    expect_true(all(is.na(logcounts(output)[is.stuff,])))

    alt <- sce
    logcounts(altExp(alt)) <- log2(counts(altExp(alt)) + 1)
    output <- unsplitAltExps(alt, prefix.rows=FALSE)
    expect_true(all(is.na(logcounts(output)[is.first,])))
    expect_identical(as.matrix(assay(altExp(alt), "logcounts")), as.matrix(logcounts(output)[is.spike,]))
    expect_true(all(is.na(logcounts(output)[is.stuff,])))

    # Counterpart is present in all assays.
    alt <- sce
    logcounts(alt) <- log2(counts(alt) + 1)
    logcounts(altExp(alt)) <- log2(counts(altExp(alt)) + 1)
    assay(altExp(alt, 2), "logcounts") <- log2(assay(altExp(alt, 2), "counts") + 1)
    output <- unsplitAltExps(alt, prefix.row=FALSE)
    expect_identical(as.matrix(logcounts(output)), rbind(logcounts(alt), assay(altExp(alt), "logcounts"), assay(altExp(alt, 2), "logcounts")))
})

test_that("rowRanges merges work as expected", {
    replacement <- rep(List(GRanges("chrA:1-100")), nrow(sce))
    names(replacement) <- rownames(sce)
    rowRanges(sce) <- replacement
    output <- unsplitAltExps(sce)

    len <- lengths(rowRanges(output))
    expect_identical(expected_names, names(len))
    expect_identical(unname(len), c(rep(1L, nrow(sce)), integer(nrow(altExp(sce))), integer(nrow(altExp(sce, 2)))))

    # What happens if one of them is a GRangesList, and another is a GRanges?
    replacement <- GRanges(rep("chrB:1-100", nrow(altExp(sce))))
    names(replacement) <- rownames(altExp(sce))
    rowRanges(altExp(sce)) <- replacement
    suppressWarnings(output <- unsplitAltExps(sce))

    expect_s4_class(rowRanges(output), "GRangesList")
    is.A <- any(seqnames(rowRanges(output))=="chrA")
    expect_true(all(is.A[is.first]))
    expect_true(!any(is.A[-is.first]))
    is.B <- any(seqnames(rowRanges(output))=="chrB")
    expect_true(all(is.B[is.spike]))
    expect_true(!any(is.B[-is.spike]))
    expect_true(all(lengths(rowRanges(output))[is.stuff]==0))

    # What happens if all of them are GRanges?
    replacement <- GRanges(rep("chrA:1-100", nrow(sce)))
    names(replacement) <- rownames(sce)
    rowRanges(sce) <- replacement

    suppressWarnings(output <- unsplitAltExps(sce))
    expect_s4_class(rowRanges(output), "GRanges")
    expect_identical(as.character(seqnames(output)), 
        rep(c("chrA", "chrB", "unknown"), c(nrow(sce), nrow(altExp(sce)), nrow(altExp(sce, 2)))))

    # How are metadata fields handled?
    rowData(sce)$thingy <- 1
    rowData(altExp(sce))$blah <- "A"
    rowData(altExp(sce, 2))$whee <- TRUE

    suppressWarnings(output <- unsplitAltExps(sce))
    expect_type(rowData(output)$thingy, "double")
    expect_type(rowData(output)$blah, "character")
    expect_type(rowData(output)$whee, "logical")
})

test_that("colData merges work as expected", {
    output <- unsplitAltExps(sce)
    expect_identical(output$Cell_Cycle, sce$Cell_Cycle)
    expect_identical(output$Protein.something, altExp(sce, 2)$something)

    output <- unsplitAltExps(sce, prefix.cols=FALSE)
    expect_identical(output$something, altExp(sce, 2)$something)
})

test_that("reducedDim merges work as expected", {
    reducedDim(sce, "PCA") <- matrix(rnorm(ncol(sce)*2), ncol=2)
    output <- unsplitAltExps(sce)
    expect_identical(reducedDims(sce), reducedDims(output))

    reducedDim(altExp(sce), "PCA") <- matrix(rnorm(ncol(sce)*2), ncol=2)
    output <- unsplitAltExps(sce)
    expect_identical(reducedDimNames(output), c("PCA", "Spike.PCA"))
    expect_identical(reducedDim(output, "Spike.PCA"), reducedDim(altExp(sce)))

    output <- unsplitAltExps(sce, prefix.cols=FALSE)
    expect_identical(reducedDimNames(output), c("PCA", "PCA"))
    expect_identical(reducedDim(output, 2), reducedDim(altExp(sce)))
})

test_that("the effect of splitAltExps is reversed", {
    alt <- sce
    altExps(alt) <- NULL
    ref <- rbind(alt, altExp(sce))
    int_metadata(ref) <- int_metadata(alt)

    alt <- sce
    altExps(alt) <- altExps(alt)[1]
    out <- unsplitAltExps(alt, delayed=FALSE, prefix.cols=NA, prefix.rows=FALSE)
    expect_identical(ref, out)
})

