#' @title
#' Miscellaneous SingleCellExperiment methods
#'
#' @description
#' Miscellaneous methods for the \linkS4class{SingleCellExperiment} class
#' that do not fit in any other documentation category.
#' 
#' @section Available methods:
#' In the following code snippets, \code{x} and \code{object} are \linkS4class{SingleCellExperiment} objects.
#' \describe{
#' \item{\code{show(object)}:}{Print a message to screen describing the contents of \code{object}.}
#' \item{\code{objectVersion(x)}:}{Return the version of the package with which \code{x} was constructed.}
#' \item{\code{sizeFactors(object)}:}{Return a numeric vector of size factors of length equal to \code{ncol(object)}.
#' If no size factors are available in \code{object}, return \code{NULL} instead.}
#' \item{\code{sizeFactors(object) <- value}:}{Replace the size factors with \code{value},
#' usually expected to be a numeric vector or vector-like object.
#' Alternatively, \code{value} can be \code{NULL} in which case any size factors in \code{object} are removed.}
#' }
#' 
#' @author Aaron Lun
#' @seealso
#' \code{\link{updateObject}}, where \code{objectVersion} is used.
#' 
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#'
#' show(sce)
#' 
#' objectVersion(sce)
#' 
#' # Setting/getting size factors.
#' sizeFactors(sce) <- runif(ncol(sce))
#' sizeFactors(sce)
#'
#' sizeFactors(sce) <- NULL
#' sizeFactors(sce)
#'
#' @name SCE-miscellaneous
#' @rdname miscellaneous
#' @docType methods
#' @aliases
#' show,SingleCellExperiment-method
#' objectVersion
#' objectVersion,SingleCellExperiment-method
NULL

#' @export
setMethod("objectVersion", "SingleCellExperiment", function(x) {
    int_metadata(x)$version
})

#' @importFrom S4Vectors coolcat
.sce_show <- function(object) {
    callNextMethod()
    coolcat("reducedDimNames(%d): %s\n", reducedDimNames(object))

    me <- mainExpName(object)
    if (is.null(me)) me <- "NULL"
    cat(sprintf("mainExpName: %s\n", me))

    coolcat("altExpNames(%d): %s\n", altExpNames(object))
}

#' @export
setMethod("show", "SingleCellExperiment", .sce_show)
