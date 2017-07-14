# Checks the combining methods.
# library(SingleCellExperiment); library(testthat); 

set.seed(1000)
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
pca <- matrix(runif(ncells*5), ncells)
sce <- SingleCellExperiment(assay=v, reducedDims=SimpleList(PCA=pca)) # A fully loaded object.
isSpike(sce, "ERCC") <- rbinom(nrow(v), 1, 0.2)==1
sizeFactors(sce) <- 2^rnorm(ncells)

# rbind
shuffled <- sample(nrow(v))
sce.alt <- sce[shuffled,]

sce2 <- rbind(sce, sce.alt)
expect_equivalent(assay(sce2), rbind(assay(sce), assay(sce.alt)))
expect_identical(isSpike(sce2), c(isSpike(sce), isSpike(sce.alt)))
expect_identical(sizeFactors(sce2), sizeFactors(sce))
expect_identical(reducedDims(sce2), reducedDims(sce))

sce.err <- sce
isSpike(sce.err, "ERCC") <- NULL
expect_error(rbind(sce.err, sce), "DataFrame 1 does not have 'is_spike_ERCC', 'is_spike'")
expect_error(rbind(sce, sce.err), "DataFrame 2 does not have 'is_spike_ERCC', 'is_spike'")

sce.lost <- sce
sizeFactors(sce.lost) <- NULL
sce2 <- rbind(sce.lost, sce)
expect_identical(sizeFactors(sce2), NULL)
sce2 <- rbind(sce, sce.lost)
expect_identical(sizeFactors(sce2), sizeFactors(sce))

sce.lost <- sce
reducedDim(sce.lost, "PCA") <- NULL
sce2 <- rbind(sce.lost, sce)
expect_equivalent(reducedDims(sce2), SimpleList())
sce2 <- rbind(sce, sce.lost)
expect_identical(reducedDims(sce2), reducedDims(sce))

# cbind
shuffled <- sample(ncells)
sce.alt <- sce[,shuffled]

sce2 <- cbind(sce, sce.alt)
expect_equivalent(assay(sce2), cbind(assay(sce), assay(sce.alt)))
expect_identical(sizeFactors(sce2), c(sizeFactors(sce), sizeFactors(sce.alt)))
expect_identical(reducedDim(sce2, "PCA"), rbind(reducedDim(sce, "PCA"), reducedDim(sce.alt, "PCA")))
expect_identical(isSpike(sce2), isSpike(sce))

sce.err <- sce
sizeFactors(sce.err) <- NULL
expect_error(cbind(sce.err, sce), "DataFrame 1 does not have 'size_factor'")
expect_error(cbind(sce, sce.err), "DataFrame 2 does not have 'size_factor'")

sce.err <- sce
reducedDim(sce.err, "PCA") <- NULL
expect_error(cbind(sce.err, sce), "object 1 does not have 'PCA' in 'reducedDims'")
expect_error(cbind(sce, sce.err), "object 2 does not have 'PCA' in 'reducedDims'")

sce.lost <- sce
isSpike(sce.lost, "ERCC") <- NULL
sce2 <- cbind(sce.lost, sce)
expect_identical(isSpike(sce2), NULL)
sce2 <- cbind(sce, sce.lost)
expect_identical(isSpike(sce2), isSpike(sce))


