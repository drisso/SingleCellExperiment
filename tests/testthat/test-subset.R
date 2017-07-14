# Checks the subsetting methods.
# library(SingleCellExperiment); library(testthat); source("test-subset.R")

set.seed(1000)
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
pca <- matrix(runif(ncells*5), ncells)
sce <- SingleCellExperiment(assay=v, reducedDims=SimpleList(PCA=pca)) # A fully loaded object.
isSpike(sce, "ERCC") <- rbinom(nrow(v), 1, 0.2)==1
sizeFactors(sce) <- 2^rnorm(ncells)

# Subsetting by row.
rownames(sce) <- paste0("Gene", seq_len(nrow(v))) 
rowData(sce)$indicator <- seq_len(nrow(sce))

for (i in 1:3) {
    if (i==1L) { 
        by.row <- sample(nrow(v), 20)
        sub.sce <- sce[by.row,]
        expect_identical(rowData(sub.sce)$indicator, by.row)
    } else if (i==2L) {
        by.row <- rbinom(nrow(v), 1, 0.2)==1
        sub.sce <- sce[by.row,]
        expect_identical(rowData(sub.sce)$indicator, which(by.row))
    } else if (i==3L) {
        by.row <- rownames(sce)[sample(nrow(v), 100)]
        sub.sce <- sce[by.row,]
        expect_identical(rowData(sub.sce)$indicator, match(by.row, rownames(sce)))
    }
    ind <- rowData(sub.sce)$indicator

    expect_identical(assay(sce)[ind,,drop=FALSE], assay(sub.sce)) # check SE elements are subsetted.
    expect_identical(isSpike(sub.sce), isSpike(sce)[ind])

    # Unchanged elements:
    expect_identical(sizeFactors(sub.sce), sizeFactors(sce)) 
    expect_identical(reducedDims(sub.sce), reducedDims(sce)) 
    expect_identical(spikeNames(sub.sce), spikeNames(sce))
    expect_identical(objectVersion(sub.sce), objectVersion(sce))
}

expect_error(sce[nrow(sce)+1,], "subscript contains out-of-bounds indices", fixed=TRUE)
expect_error(sce["A",], "index out of bounds: A")

# Subsetting by column.
colnames(sce) <- paste0("Cell", seq_len(ncells))
colData(sce)$indicator <- seq_len(ncells)

for (j in 1:3) {
    if (j==1L) { 
        by.col <- sample(ncells, 20)
        sub.sce <- sce[,by.col]
        expect_identical(colData(sub.sce)$indicator, by.col)
    } else if (j==2L) {
        by.col <- rbinom(ncells, 1, 0.2)==1
        sub.sce <- sce[,by.col]
        expect_identical(colData(sub.sce)$indicator, which(by.col))
    } else if (j==3L) {
        by.col <- colnames(sce)[sample(ncells, 50)]
        sub.sce <- sce[,by.col]
        expect_identical(colData(sub.sce)$indicator, match(by.col, colnames(sce)))
    }
    ind <- colData(sub.sce)$indicator

    expect_identical(assay(sce)[,ind,drop=FALSE], assay(sub.sce)) # check SE elements are subsetted.
    expect_identical(sizeFactors(sub.sce), sizeFactors(sce)[ind])
    expect_identical(reducedDim(sub.sce, "PCA"), pca[ind,,drop=FALSE])
 
    # Unchanged elements:
    expect_identical(isSpike(sub.sce), isSpike(sce))
    expect_identical(isSpike(sub.sce), isSpike(sce))
    expect_identical(spikeNames(sub.sce), spikeNames(sce))
    expect_identical(objectVersion(sub.sce), objectVersion(sce))
}

expect_error(sce[,ncells+1], "subscript contains out-of-bounds indices", fixed=TRUE)
expect_error(sce[,"A"], "index out of bounds: A")

# Subset replacement by row.
sce.alt <- sce
rownames(sce.alt) <- paste0(rownames(sce), "x")
SingleCellExperiment:::int_metadata(sce.alt)$whee <- 1

scex <- sce.alt
to <- 1:10
from <- 21:30
scex[to,] <- sce[from,]

expect_identical(assay(scex)[to,,drop=FALSE], assay(sce)[from,,drop=FALSE])
expect_equivalent(assay(scex)[-to,,drop=FALSE], assay(sce)[-to,,drop=FALSE])
expect_identical(isSpike(scex)[to], isSpike(sce)[from])
expect_identical(isSpike(scex)[-to], isSpike(sce)[-to])
expect_identical(SingleCellExperiment:::int_metadata(scex), SingleCellExperiment:::int_metadata(sce))
expect_identical(sizeFactors(scex), sizeFactors(sce))

scex2 <- sce.alt
isSpike(scex2, "ERCC") <- NULL
expect_error(scex2[to,] <- sce[from,], "DataFrame 1 does not have 'is_spike_ERCC', 'is_spike'")
scex2 <- scex3 <- sce.alt
isSpike(scex3, "ERCC") <- NULL
expect_error(scex2[to,] <- scex3[from,], "DataFrame 2 does not have 'is_spike_ERCC', 'is_spike'")

scex2 <- sce.alt # Again for character
to <- rownames(scex2)[1:10]
from <- rownames(sce)[21:30]
scex2[to,] <- sce[from,]
expect_equal(scex, scex2)

# Subset replacement by column.
sce.alt <- sce
colnames(sce.alt) <- paste0(colnames(sce), "x")

scex <- sce.alt
to <- 1:10
from <- 21:30
scex[,to] <- sce[,from]

expect_identical(assay(scex)[,to,drop=FALSE], assay(sce)[,from,drop=FALSE])
expect_equivalent(assay(scex)[,-to,drop=FALSE], assay(sce)[,-to,drop=FALSE])
expect_identical(sizeFactors(scex)[to], sizeFactors(sce)[from])
expect_identical(sizeFactors(scex)[-to], sizeFactors(sce)[-to])
expect_identical(SingleCellExperiment:::int_metadata(scex), SingleCellExperiment:::int_metadata(sce))
expect_identical(isSpike(scex), isSpike(sce))

scex2 <- sce.alt
sizeFactors(scex2) <- NULL
expect_error(scex2[,to] <- sce[,from], "DataFrame 1 does not have 'size_factor'")
scex2 <- scex3 <- sce.alt
sizeFactors(scex3) <- NULL
expect_error(scex2[,to] <- scex3[,from], "DataFrame 2 does not have 'size_factor'")
scex2 <- sce.alt
reducedDim(scex2, "PCA") <- NULL
expect_error(scex2[,to] <- sce[,from], "object 1 does not have 'PCA' in 'reducedDims'")

scex2 <- sce.alt # Again for character
to <- colnames(scex2)[1:10]
from <- colnames(sce)[21:30]
scex2[,to] <- sce[,from]
expect_equal(scex, scex2)

# Subset replacement, all.
sce.alt <- sce
rownames(sce.alt) <- paste0(rownames(sce), "x")
colnames(sce.alt) <- paste0(colnames(sce), "x")
SingleCellExperiment:::int_metadata(sce.alt)$whee <- 1

scex <- sce.alt
scex[] <- sce
expect_equal(scex, sce)


