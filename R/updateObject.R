#' @export
#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors DataFrame
setMethod("updateObject", "SingleCellExperiment", function(object, ..., verbose=FALSE) {
    if (objectVersion(object) < "1.7.1") {
        stuff <- object@reducedDims                        
        int_colData(object)$reducedDims <- do.call(DataFrame, lapply(stuff, I))

        if (verbose) {
            message("[updateObject] ", class(object), " object uses ", 
                "internal representation\n", "[updateObject] from SingleCellExperiment ", 
                objectVersion(object), ". ", "Updating it ... ", appendLF = FALSE)
        }
    }
    object
})
