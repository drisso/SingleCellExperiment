#' @export
#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors DataFrame
setMethod("updateObject", "SingleCellExperiment", function(object, ..., verbose=FALSE) {
    if (objectVersion(object) < "1.7.1") {
        stuff <- object@reducedDims
        reducedDims(object) <- stuff
        altExps(object) <- NULL

        if (verbose) {
            message("[updateObject] ", class(object), " object uses ", 
                "internal representation\n", "[updateObject] from SingleCellExperiment ", 
                objectVersion(object), ". ", "Updating it ... ", appendLF = FALSE)
        }
    }
    object
})
