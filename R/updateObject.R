#' Update a SingleCellExperiment object
#'
#' Update \linkS4class{SingleCellExperiment} objects to the latest version of the class structure.
#' This is usually called by methods in the \pkg{SingleCellExperiment} package rather than by users or downstream packages.
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
#' \item Objects created before 1.9.1 are modified so that the size factors are stored by \code{\link{sizeFactors<-}} in \code{\link{colData}} rather than \code{\link{int_colData}}.
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
    old.ver <- objectVersion(object)
    triggered <- FALSE

    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    # Update possibly outdated DataFrame instances.
    if (.hasSlot(object, "int_elementMetadata"))
        object@int_elementMetadata <-
            updateObject(object@int_elementMetadata, ..., verbose=verbose)
    if (.hasSlot(object, "int_colData"))
        object@int_colData <-
            updateObject(object@int_colData, ..., verbose=verbose)

    if (old.ver < "1.7.1") {
        # Need this to avoid a circular recursion when calling reducedDims()<-.
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

        triggered <- TRUE
    }

    if (old.ver < "1.9.2") {
        # Need this to avoid a circular recursion when calling sizeFactors()<-.
        int_metadata(object)$version <- packageVersion("SingleCellExperiment")

        old_sf <- int_colData(object)$size_factor
        has_sf <- sizeFactors(object)
        if (!is.null(has_sf) && !identical(has_sf, old_sf)) {
            warning(sprintf("clobbering old 'sizeFactors' in 'colData(<%s>)'", class(object)[1]))
        }
        sizeFactors(object) <- old_sf

        triggered <- TRUE
    }

    if (old.ver < "1.11.3") {
        # Need this to avoid a circular recursion when calling rowPairs()<-.
        int_metadata(object)$version <- packageVersion("SingleCellExperiment")
        
        rowPairs(object) <- list()
        colPairs(object) <- list()

        triggered <- TRUE
    }

    if (verbose && triggered) {
        message("[updateObject] ", class(object)[1], " object uses ", 
            "internal representation\n", "[updateObject] from SingleCellExperiment ", 
            old.ver, ". ", "Updating it ...\n", appendLF = FALSE)
    }

    int_metadata(object)$version <- packageVersion("SingleCellExperiment")
    callNextMethod()
})
