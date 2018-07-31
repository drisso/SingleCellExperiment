# Internal functions not intended for users.

setMethod("int_elementMetadata", "SingleCellExperiment", function(x) x@int_elementMetadata)
setReplaceMethod("int_elementMetadata", "SingleCellExperiment", function(x, value) {
    x@int_elementMetadata <- value
    return(x)
})

setMethod("int_colData", "SingleCellExperiment", function(x) x@int_colData)
setReplaceMethod("int_colData", "SingleCellExperiment", function(x, value) {
    x@int_colData <- value
    return(x)
})

setMethod("int_metadata", "SingleCellExperiment", function(x) x@int_metadata)
setReplaceMethod("int_metadata", "SingleCellExperiment", function(x, value) {
    x@int_metadata <- value
    return(x)
})

#' @export
setMethod("objectVersion", "SingleCellExperiment", function(x) {
    int_metadata(x)$version
})
