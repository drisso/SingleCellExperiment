# Checks for proper functioning of the methods.
# library(SingleCellExperiment); library(testthat); source("test-methods.R")

set.seed(1000)
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
sce <- SingleCellExperiment(assay=v)

# Adding spike-ins.

test_that("spike-in getters/setters are functioning", {
    is.spike1 <- rbinom(nrow(v), 1, 0.2)==1
    isSpike(sce, "ERCC") <- is.spike1
    expect_identical(spikeNames(sce), "ERCC")
    expect_identical(isSpike(sce, "ERCC"), is.spike1)
    expect_identical(isSpike(sce), is.spike1)
    
    is.spike2 <- rbinom(nrow(v), 1, 0.3)==1
    isSpike(sce, "SIRV") <- is.spike2
    expect_identical(spikeNames(sce), c("ERCC", "SIRV"))
    expect_identical(isSpike(sce, "ERCC"), is.spike1) # check still the same.
    expect_identical(isSpike(sce, "SIRV"), is.spike2)
    expect_identical(isSpike(sce), is.spike1 | is.spike2)
    
    isSpike(sce, "ERCC") <- NULL
    expect_identical(spikeNames(sce), "SIRV")
    expect_error(isSpike(sce, "ERCC"), "spike-in set 'ERCC' does not exist")
    expect_identical(isSpike(sce, "SIRV"), is.spike2)
    expect_identical(isSpike(sce), is.spike2)
    
    chosen <- sample(nrow(v), 20) # integer setting
    isSpike(sce, "ERCC") <- chosen
    expect_identical(which(isSpike(sce, "ERCC")), sort(chosen))
    expect_identical(spikeNames(sce), c("SIRV", "ERCC")) # flipped
    
    rownames(sce) <- paste0("Gene", seq_len(nrow(v))) # character setting
    isSpike(sce, "SIRV") <- rownames(sce)[chosen]
    expect_identical(which(isSpike(sce, "SIRV")), sort(chosen))
    rownames(sce) <- NULL
    
    SingleCellExperiment:::int_metadata(sce)$spike_names <- c("random")
    expect_error(validObject(sce), "no field specifying rows belonging to spike-in set 'random'", fixed=TRUE)
})

# Adding size factors.

test_that("size factor getters/setters are functioning", {
    sf1 <- 2^rnorm(ncells)
    sizeFactors(sce) <- sf1
    expect_identical(sizeFactors(sce), sf1)
    
    sf2 <- 2^rnorm(ncells, sd=2)
    sizeFactors(sce, "ERCC") <- sf2
    expect_identical(sizeFactors(sce), sf1) # check still the same
    expect_identical(sizeFactors(sce, "ERCC"), sf2)
    
    sizeFactors(sce) <- NULL
    expect_identical(sizeFactors(sce), NULL)
    expect_identical(sizeFactors(sce, "ERCC"), sf2) # check still the same
    
    sizeFactors(sce, "ERCC") <- NULL
    expect_identical(sizeFactors(sce, "ERCC"), NULL)
})

# Adding reduced dimensions.

test_that("reduced dimension getters/setters are functioning", {
    d1 <- matrix(rnorm(ncells*4), ncol=4)
    d2 <- matrix(rnorm(ncells), ncol=1)

    reducedDim(sce, "PCA") <- d1
    expect_identical(reducedDim(sce, "PCA"), d1)
    expect_identical(reducedDims(sce), SimpleList(PCA=d1))
    expect_identical(reducedDimNames(sce), "PCA")

    reducedDim(sce, "tSNE") <- d2
    expect_identical(reducedDim(sce, "tSNE"), d2)
    expect_identical(reducedDims(sce), SimpleList(PCA=d1, tSNE=d2))
    expect_identical(reducedDimNames(sce), c("PCA", "tSNE"))

    reducedDim(sce, "PCA") <- NULL
    expect_identical(reducedDim(sce, "PCA"), NULL)
    expect_identical(reducedDim(sce, "tSNE"), d2)
    expect_identical(reducedDims(sce), SimpleList(tSNE=d2))
    expect_identical(reducedDimNames(sce), "tSNE")

    reducedDims(sce) <- SimpleList()
    expect_identical(reducedDims(sce), SimpleList())
    expect_identical(reducedDim(sce), NULL)

    reducedDims(sce) <- SimpleList(DM=d1)
    expect_identical(SimpleList(DM=d1), reducedDims(sce))
    expect_identical(d1, reducedDim(sce))

    expect_error(reducedDim(sce, "DM") <- d1[1:10,], "each element of 'reducedDims' must be a matrix with nrow equal to 'ncol(object)'", fixed=TRUE)
    expect_error(reducedDim(sce, "DM") <- "huh", "each element of 'reducedDims' must be a matrix with nrow equal to 'ncol(object)'", fixed=TRUE)
    expect_error(reducedDim(sce, 5), "subscript is out of bounds")

    # Repeating with numeric indices, and default extraction/assignment.
    # This gets a bit confusing as the order changes when earlier elements are wiped out.

    reducedDims(sce) <- SimpleList() # Wiping.
    reducedDim(sce, 1) <- d1
    expect_identical(reducedDim(sce), d1)
    expect_identical(reducedDim(sce, 1), d1)
    reducedDim(sce, 2) <- d2
    expect_identical(reducedDim(sce), d1)
    expect_identical(reducedDim(sce, 2), d2)
    
    reducedDim(sce, 1) <- NULL # d2 becomes the first element now.
    expect_identical(reducedDim(sce), d2)
    expect_identical(reducedDim(sce, 1), d2)
    reducedDim(sce, "PCA") <- d1 # d1 is the second element.
    expect_identical(reducedDim(sce, 1), d2)
    expect_identical(reducedDim(sce, 2), d1)

    reducedDim(sce) <- NULL # d1 now becomes the first element again.
    expect_identical(reducedDim(sce), d1)
    expect_identical(reducedDimNames(sce), "PCA")
    reducedDim(sce) <- d2 # d2 now overwrites the first element.
    expect_identical(reducedDim(sce, 1), d2)
    expect_identical(reducedDimNames(sce), "PCA")

    expect_error(reducedDim(sce, 5) <- d1, "subscript is out of bounds")
})

# Checking package version.

test_that("object version extraction works", {
    expect_identical(objectVersion(sce), packageVersion("SingleCellExperiment"))
})

# Checking colData and rowData

test_that("special colData/rowData getters/setters work", {
    isSpike(sce, "ERCC") <- rbinom(nrow(v), 1, 0.2)==1
    sizeFactors(sce, "SF") <- 2^rnorm(ncells)
    sizeFactors(sce, "ERCC") <- 2^rnorm(ncells)
    
    random_coldata <- DataFrame(a=rnorm(ncells), b=runif(ncells, 0, 1))
    colData(sce) <- random_coldata
    expect_identical(colData(sce), random_coldata)
    expect_identical(colData(sce), colData(sce, internal=FALSE))
    expect_identical(colData(sce, internal=TRUE),
                     cbind(colData(sce), SingleCellExperiment:::int_colData(sce)))
    
    random_rowdata <- DataFrame(a=rnorm(NROW(sce)), b=runif(NROW(sce), 0, 1))
    rowData(sce) <- random_rowdata
    expect_identical(rowData(sce), random_rowdata)
    expect_identical(rowData(sce), rowData(sce, internal=FALSE))
    expect_identical(rowData(sce, internal=TRUE),
                     cbind(rowData(sce), SingleCellExperiment:::int_elementMetadata(sce)))
    
    rowData(sce)$is_spike <- rnorm(NROW(sce))
    expect_warning(rowData(sce, internal=TRUE), "overlapping names in internal and external rowData")
    
    colData(sce)$size_factor_ERCC <- rnorm(ncells)
    expect_warning(colData(sce, internal=TRUE), "overlapping names in internal and external colData")
})

#  Checking the assay convenience wrappers.

test_that("assay getters/setters work", {
    v2 <- matrix(runif(20000), ncol=ncells)
    counts(sce) <- v2
    expect_equivalent(counts(sce), v2)

    v3 <- log2(v2)
    logcounts(sce) <- v3
    expect_equivalent(counts(sce), v2)
    expect_equivalent(logcounts(sce), v3)

    counts(sce) <- NULL
    expect_equivalent(logcounts(sce), v3) 
    expect_error(counts(sce), "invalid subscript") 
})
