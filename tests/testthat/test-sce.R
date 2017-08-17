# Checks for proper construction and get/setting of the slots of a SingleCellExperiment.
# library(SingleCellExperiment); library(testthat)

set.seed(1000)
ncells <- 100

# Adding assays.
v <- matrix(rnorm(20000), ncol=ncells)
sce <- SingleCellExperiment(assay=v)
expect_equivalent(assay(sce), v)

u <- matrix(rpois(20000, 5), ncol=ncells)
sce <- SingleCellExperiment(assay=list(counts=u, exprs=v))
expect_equivalent(assay(sce, "counts"), u)
expect_equivalent(assay(sce, "exprs"), v)

w <- matrix(runif(20000), ncol=ncells)
assay(sce, "exprs") <- w
expect_equivalent(assay(sce, "exprs"), w)

# Adding colData and rowData
rd <- DataFrame(stuff=runif(nrow(v)))
cd <- DataFrame(whee=runif(ncells))
sce <- SingleCellExperiment(u, rowData=rd, colData=cd)
expect_equal(rd, rowData(sce))
expect_equal(cd, colData(sce))

# Coercion from SummarizedExperiment
se <- SummarizedExperiment(u, rowData = rd, colData = cd)
sce2 <- as(se, "SingleCellExperiment")
expect_equal(sce, sce2)

# Coercion from RangedSummarizedExperiment
ranges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                  IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                  strand=sample(c("+", "-"), 200, TRUE))
mcols(ranges) <- rd

rse <- SummarizedExperiment(u, colData = cd, rowRanges = ranges)
sce3 <- as(rse, "SingleCellExperiment")
expect_equal(assays(sce), assays(sce3))
expect_equal(rowData(sce), rowData(sce3))
expect_equal(colData(sce), colData(sce3))

# Manipulating metadata
cextra <- rnorm(ncells)
sce$blah <- cextra
expect_equal(cextra, colData(sce)$blah)
rextra <- rnorm(nrow(v))
rowData(sce)$blah <- rextra
expect_equal(rextra, rowData(sce)$blah)

sce <- SingleCellExperiment(u, metadata=list(yay=1))
expect_equal(metadata(sce)$yay, 1)
metadata(sce)$yay <- "stuff"
expect_identical(metadata(sce)$yay, "stuff")

# Adding reduced dimensions.
pca <- matrix(runif(ncells*5), ncells)
tsne <- matrix(rnorm(ncells*2), ncells)
combined <- SimpleList(PCA=pca, tSNE=tsne)
sce <- SingleCellExperiment(u, reducedDims=combined)

expect_equivalent(pca, reducedDim(sce, "PCA"))
expect_equivalent(tsne, reducedDim(sce, "tSNE"))
expect_equivalent(reducedDims(sce), combined)

# Checking internals.
sce <- SingleCellExperiment(assay=u)
expect_identical(nrow(SingleCellExperiment:::int_elementMetadata(sce)), nrow(sce))
expect_identical(nrow(SingleCellExperiment:::int_colData(sce)), ncol(sce))
expect_identical(length(SingleCellExperiment:::int_metadata(sce)), 1L)

SingleCellExperiment:::int_elementMetadata(sce)$whee <- rextra
expect_equal(rextra, SingleCellExperiment:::int_elementMetadata(sce)$whee)
SingleCellExperiment:::int_elementMetadata(sce) <- DataFrame(1:5)
expect_error(validObject(sce), "'nrow' of internal 'rowData' not equal to 'nrow(object)'", fixed=TRUE)

SingleCellExperiment:::int_colData(sce)$stuff <- cextra
expect_equal(cextra, SingleCellExperiment:::int_colData(sce)$stuff)
SingleCellExperiment:::int_colData(sce) <- DataFrame(1:5)
expect_error(validObject(sce), "'nrow' of internal 'colData' not equal to 'ncol(object)'", fixed=TRUE)

SingleCellExperiment:::int_metadata(sce)$urg <- "I was here"
expect_identical(SingleCellExperiment:::int_metadata(sce)$urg, "I was here")

