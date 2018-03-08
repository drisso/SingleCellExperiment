# Checks for proper functioning of the methods.
# library(SingleCellExperiment); library(testthat); source("test-methods.R")

set.seed(1000)
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
sce <- SingleCellExperiment(assay=v)

test_that("spike-in getters/setters are functioning", {
    is.spike1 <- rbinom(nrow(v), 1, 0.2)==1
    isSpike(sce, "ERCC") <- is.spike1
    expect_identical(spikeNames(sce), "ERCC")
    expect_identical(isSpike(sce, "ERCC"), is.spike1)
    expect_identical(isSpike(sce), is.spike1)

    # Check what happens when we add another spike-in set.    
    is.spike2 <- rbinom(nrow(v), 1, 0.3)==1
    isSpike(sce, "SIRV") <- is.spike2
    expect_identical(spikeNames(sce), c("ERCC", "SIRV"))
    expect_identical(isSpike(sce, "ERCC"), is.spike1) # check still the same.
    expect_identical(isSpike(sce, "SIRV"), is.spike2)
    expect_identical(isSpike(sce), is.spike1 | is.spike2)
    
    # Check what happens when we clear a spike-in set.
    isSpike(sce, "ERCC") <- NULL
    expect_identical(spikeNames(sce), "SIRV")
    expect_identical(isSpike(sce, "ERCC"), NULL)
    expect_identical(isSpike(sce, "SIRV"), is.spike2)
    expect_identical(isSpike(sce), is.spike2)
    
    # Checking that it still behaves with integers.
    chosen <- sample(nrow(v), 20) 
    isSpike(sce, "ERCC") <- chosen
    expect_identical(which(isSpike(sce, "ERCC")), sort(chosen))
    expect_identical(spikeNames(sce), c("SIRV", "ERCC")) # flipped
   
    # Checking that it behaves with character strings.
    rownames(sce) <- paste0("Gene", seq_len(nrow(v))) 
    isSpike(sce, "SIRV") <- rownames(sce)[chosen]
    expect_identical(which(isSpike(sce, "SIRV")), sort(chosen))
   
    # Checking that it throws properly if spike-in sets and spike-in fields are not in sync. 
    alt.sce <- sce 
    SingleCellExperiment:::int_metadata(alt.sce)$spike_names <- c("random")
    expect_error(validObject(alt.sce), "no field specifying rows belonging to spike-in set 'random'", fixed=TRUE)

    # Checking that clearing the spike-ins works properly.
    alt.sce <- clearSpikes(sce)
    expect_identical(spikeNames(alt.sce), character(0))
    expect_identical(isSpike(alt.sce), NULL)
    expect_identical(isSpike(alt.sce, "ERCC"), NULL)
    expect_identical(isSpike(alt.sce, "SIRV"), NULL)

    # Check that unnamed spike-ins are correctly handled.
    expect_warning(isSpike(sce) <- 1:10, "deprecated")
    expect_identical(which(isSpike(sce, "")), 1:10)
    expect_identical(spikeNames(sce), "") # for the time being.
    expect_warning(isSpike(sce) <- NULL, "deprecated")
    expect_identical(isSpike(sce), NULL)
    expect_identical(spikeNames(sce), character(0))
})

test_that("size factor getters/setters are functioning", {
    sf1 <- 2^rnorm(ncells)
    sizeFactors(sce) <- sf1
    expect_identical(sizeFactors(sce), sf1)
    expect_identical(sizeFactorNames(sce), character(0))
    
    sf2 <- 2^rnorm(ncells, sd=2)
    sizeFactors(sce, "ERCC") <- sf2
    expect_identical(sizeFactors(sce), sf1) # check still the same
    expect_identical(sizeFactors(sce, "ERCC"), sf2)
    expect_identical(sizeFactorNames(sce), "ERCC")
    
    # Automated deletion.    
    alt.sce <- clearSizeFactors(sce)
    expect_identical(sizeFactors(alt.sce), NULL)
    expect_identical(sizeFactors(alt.sce, "ERCC"), NULL)
    expect_identical(sizeFactorNames(alt.sce), character(0))

    # Checking that it throws properly if spike-in sets and spike-in fields are not in sync. 
    alt.sce <- sce 
    SingleCellExperiment:::int_metadata(alt.sce)$size_factor_names <- c("random")
    expect_error(validObject(alt.sce), "no field specifying size factors for set 'random'", fixed=TRUE)

    # Manual deletion.
    sizeFactors(sce) <- NULL
    expect_identical(sizeFactors(sce), NULL)
    expect_identical(sizeFactors(sce, "ERCC"), sf2) # check still the same
    
    sizeFactors(sce, "ERCC") <- NULL
    expect_identical(sizeFactors(sce, "ERCC"), NULL)
})

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

    expect_error(reducedDim(sce, 5) <- d1, "subscript is out of bounds")
})

test_that("object version extraction works", {
    expect_identical(objectVersion(sce), packageVersion("SingleCellExperiment"))
})

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

test_that("assay getters/setters work", {
    v2 <- matrix(runif(20000), ncol=ncells)
    counts(sce) <- v2
    expect_equivalent(counts(sce), v2)

    v3 <- log2(v2)
    logcounts(sce) <- v3
    expect_equivalent(counts(sce), v2)
    expect_equivalent(logcounts(sce), v3)

    cpm(sce) <- v3 + v2
    expect_equivalent(cpm(sce), v3+v2)
    tpm(sce) <- v3 - v2
    expect_equivalent(tpm(sce), v3-v2)

    counts(sce) <- NULL
    expect_equivalent(logcounts(sce), v3) 
    expect_error(counts(sce), "invalid subscript") 
})
