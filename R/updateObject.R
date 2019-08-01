#' @export
#' @importFrom BiocGenerics updateObject
#' @importFrom utils packageVersion
setMethod("updateObject", "SingleCellExperiment", function(object, ..., verbose=FALSE) {
    if (objectVersion(object) < "1.7.1") {
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }

        int_metadata(object)$version <- packageVersion("SingleCellExperiment")
        reducedDims(object) <- object@reducedDims
        altExps(object) <- NULL

        if (verbose) {
            message("[updateObject] ", class(object), " object uses ", 
                "internal representation\n", "[updateObject] from SingleCellExperiment ", 
                objectVersion(object), ". ", "Updating it ... ", appendLF = FALSE)
        }
    }
    object
})
