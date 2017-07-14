# Checks for proper functioning of the methods.
# library(SingleCellExperiment); library(testthat)

set.seed(1000)
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
sce <- SingleCellExperiment(assay=v)

# Adding spike-ins.
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

# Adding size factors.
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

# Checking package version.
expect_identical(objectVersion(sce), packageVersion("SingleCellExperiment"))

# Checking colData and rowData
sizeFactors(sce, "SF") <- sf1
sizeFactors(sce, "ERCC") <- sf2

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
expect_warning(rowData(sce, internal=TRUE), "Overlapping column names")

colData(sce)$size_factor_ERCC <- rnorm(ncells)
expect_warning(colData(sce, internal=TRUE), "Overlapping column names")
