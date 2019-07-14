# Basic methods for the SummarizedExperimentByColumn class,
# an internal class that powers the alt* suite of methods.

.get_se <- function(x) x@se

SummarizedExperimentByColumn <- function(se) new("SummarizedExperimentByColumn", se=se)

setMethod("length", "SummarizedExperimentByColumn", function(x) ncol(.get_se(x)))

setMethod("[", "SummarizedExperimentByColumn", function(x, i, j, ..., drop=FALSE) {
    initialize(x, se=.get_se(x)[,i])
})

setReplaceMethod("[", "SummarizedExperimentByColumn", function(x, i, j, ..., value) {
    left <- .get_se(x)
    left[,i] <- .get_se(value)
    initialize(x, se=left)
})

setMethod("c", "SummarizedExperimentByColumn", function(x, ...) {
    gathered <- lapply(list(x, ...), .get_se)
    initialize(x, se=do.call(cbind, gathered))
})

setMethod("names", "SummarizedExperimentByColumn", function(x) colnames(.get_se(x)))

setReplaceMethod("names", "SummarizedExperimentByColumn", function(x, value) {
    colnames(x@se) <- value
    x
})
