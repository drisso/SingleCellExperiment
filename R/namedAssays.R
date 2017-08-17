# This defines some convenience wrappers for common entires in the assays slot.

GET_FUN <- function(exprs_values) {
    (exprs_values) # To ensure evaluation
    function(object) {
        assay(object, i=exprs_values)
    }
}

SET_FUN <- function(exprs_values) {
    (exprs_values) # To ensure evaluation
    function(object, value) {
        assay(object, i=exprs_values) <- value
        object
    }
}

for (x in c("counts", "normcounts", "logcounts", "cpm", "tpm")) { 
    setMethod(x, "SingleCellExperiment", GET_FUN(x))
    setReplaceMethod(x, c("SingleCellExperiment", "ANY"), SET_FUN(x))
}
