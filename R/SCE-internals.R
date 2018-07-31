# Internal functions not intended for users.

#' @export
setMethod("int_elementMetadata", "SingleCellExperiment", function(x) x@int_elementMetadata)

#' @export
setReplaceMethod("int_elementMetadata", "SingleCellExperiment", function(x, value) {
    x@int_elementMetadata <- value
    return(x)
})

#' @export
setMethod("int_colData", "SingleCellExperiment", function(x) x@int_colData)

#' @export
setReplaceMethod("int_colData", "SingleCellExperiment", function(x, value) {
    x@int_colData <- value
    return(x)
})

#' @export
setMethod("int_metadata", "SingleCellExperiment", function(x) x@int_metadata)

#' @export
setReplaceMethod("int_metadata", "SingleCellExperiment", function(x, value) {
    x@int_metadata <- value
    return(x)
})

