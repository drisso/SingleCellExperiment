# Setting up the options for a fully loaded SingleCellExperiment.

set.seed(1000)
ncells <- 100
u <- matrix(rpois(20000, 5), ncol=ncells)
v <- matrix(rnorm(20000), ncol=ncells)
empty <- SingleCellExperiment(list(counts=u, logcounts=v))
loaded <- empty

#########################################
# Mocking up some reduced dimensions.

d1 <- matrix(rnorm(ncells*4), ncol=4)
d2 <- matrix(rnorm(ncells*2), ncol=2)
reducedDim(loaded, "PCA") <- d1
reducedDim(loaded, "TSNE") <- d2

#########################################
# Mocking up some alternative Experiments.

se1 <- SummarizedExperiment(
    list(counts=matrix(rpois(1000, 5), ncol=ncells),
        logcounts=matrix(rnorm(1000, 5), ncol=ncells)
    )
)
rowData(se1)$stuff <- sample(LETTERS, nrow(se1), replace=TRUE)
rownames(se1) <- sprintf("SPIKE_%i", seq_len(nrow(se1)))

se2 <- SummarizedExperiment(
    list(counts=matrix(rpois(500, 5), ncol=ncells),
        logcounts=matrix(rnorm(500), ncol=ncells)
    )
)
rowData(se2)$blah <- sample(letters, nrow(se2), replace=TRUE)
rownames(se2) <- sprintf("TAG_%i", seq_len(nrow(se2)))

altExperiment(loaded, "Spike") <- se1
altExperiment(loaded, "Protein") <- se2

#########################################
# Other load-ups. 

isSpike(loaded, "ERCC") <- rbinom(nrow(v), 1, 0.2)==1
sizeFactors(loaded) <- 2^rnorm(ncells)
