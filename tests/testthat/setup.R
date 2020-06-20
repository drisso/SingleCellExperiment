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

altExp(loaded, "Spike") <- se1
altExp(loaded, "Protein") <- se2

#########################################
# Adding some pairings.

rhits <- SelfHits(
    sample(nrow(loaded), 100),
    sample(nrow(loaded), 100),
    nnode=nrow(loaded)
)
mcols(rhits)$value <- runif(100)

chits <- SelfHits(
    sample(ncol(loaded), 20),
    sample(ncol(loaded), 20),
    nnode=ncol(loaded)
)
mcols(chits)$value <- runif(20)

rowPair(loaded) <- rhits
colPair(loaded) <- chits

#########################################
# Other load-ups. 

sizeFactors(loaded) <- 2^rnorm(ncells)
