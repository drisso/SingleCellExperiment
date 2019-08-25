#' Update a SingleCellExperiment object
#' 
#' @param object A old \linkS4class{SingleCellExperiment} object.
#' @param ... Additional arguments that are ignored.
#' @param verbose Logical scalar indicating whether a message should be emitted as the object is updated.
#' 
#' @details
#' This function updates the SingleCellExperiment to match changes in the internal class representation.
#' Changes are as follows:
#' \itemize{
#' \item Objects created before 1.7.1 are modified to include \code{\link{altExps}} and \code{\link{reducedDims}} fields in their internal column metadata.
#' Reduced dimension results previously in the \code{reducedDims} slot are transferred to the \code{reducedDims} field.
#' }
#' 
#' @return
#' An updated version of \code{object}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{objectVersion}}, which is used to determine if the object is up-to-date.
#' 
#' @name updateObject 
#' @export
#' @aliases updateObject updateObject,SingleCellExperiment-method
#' @importFrom BiocGenerics updateObject
#' @importFrom S4Vectors DataFrame
#' @importFrom utils packageVersion
setMethod("updateObject", "SingleCellExperiment", function(object, ..., verbose=FALSE) {
    if (objectVersion(object) < "1.7.1") {
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }

        int_metadata(object)$version <- packageVersion("SingleCellExperiment")

        if (!.red_key %in% colnames(int_colData(object))) {
            reddims <- NULL

            # Need the try() to handle rare cases involving an as() from old
            # subclass, which does not have @reducedDims anymore. See
            # https://github.com/drisso/SingleCellExperiment/issues/37.
            try(reddims <- object@reducedDims, silent=TRUE) 
            reducedDims(object) <- reddims
        }

        if (!.alt_key %in% colnames(int_colData(object))) { 
            altExps(object) <- NULL
        }

        if (verbose) {
            message("[updateObject] ", class(object), " object uses ", 
                "internal representation\n", "[updateObject] from SingleCellExperiment ", 
                objectVersion(object), ". ", "Updating it ... ", appendLF = FALSE)
        }
    }

    int_metadata(object)$version <- packageVersion("SingleCellExperiment")
    object
})
