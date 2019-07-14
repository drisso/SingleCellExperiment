#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors DataFrame
setMethod("updateObject", "SingleCellExperiment", function(object, ..., verbose=FALSE) {
    if (objectVersion(object) < "1.5.3") {
        stuff <- object@reducedDims                        
        int_colData(object)$reducedDims <- do.call(DataFrame, lapply(stuff, I))
    }
    object
})
