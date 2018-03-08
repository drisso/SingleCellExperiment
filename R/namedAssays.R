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

#' @export
#' @importFrom BiocGenerics counts
setMethod("counts", "SingleCellExperiment", GET_FUN("counts"))

#' @export
#' @importFrom BiocGenerics "counts<-"
setReplaceMethod("counts", c("SingleCellExperiment", "ANY"), SET_FUN("counts"))

#' @export
setMethod("logcounts", "SingleCellExperiment", GET_FUN("logcounts"))

#' @export
setReplaceMethod("logcounts", c("SingleCellExperiment", "ANY"), SET_FUN("logcounts"))

#' @export
setMethod("normcounts", "SingleCellExperiment", GET_FUN("normcounts"))

#' @export
setReplaceMethod("normcounts", c("SingleCellExperiment", "ANY"), SET_FUN("normcounts"))

#' @export
setMethod("cpm", "SingleCellExperiment", GET_FUN("cpm"))

#' @export
setReplaceMethod("cpm", c("SingleCellExperiment", "ANY"), SET_FUN("cpm"))

#' @export
setMethod("tpm", "SingleCellExperiment", GET_FUN("tpm"))

#' @export
setReplaceMethod("tpm", c("SingleCellExperiment", "ANY"), SET_FUN("tpm"))
