#' @export
#' @importFrom BiocGenerics updateObject
setMethod("updateObject", "SingleCellExperiment", function(object, ..., verbose=FALSE) {
    if (objectVersion(object) < "1.7.1") {
        stuff <- object@reducedDims
        object <- .rse_to_sce(as(object, "RangedSummarizedExperiment"), reducedDims=stuff)

        if (verbose) {
            message("[updateObject] ", class(object), " object uses ", 
                "internal representation\n", "[updateObject] from SingleCellExperiment ", 
                objectVersion(object), ". ", "Updating it ... ", appendLF = FALSE)
        }
    }
    object
})
