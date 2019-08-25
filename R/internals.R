#' Internal SingleCellExperiment functions 
#' 
#' Methods to get or set internal fields from the SingleCellExperiment class.
#' Thse functions are intended for package developers who want to add protected fields to a SingleCellExperiment.
#' They should \emph{not} be used by ordinary users of the \pkg{SingleCellExperiment} package.
#' 
#' @section Getters:
#' In the following code snippets, \code{x} is a \linkS4class{SingleCellExperiment}.
#' \describe{
#' \item{\code{int_elementMetadata(x)}:}{Returns a \linkS4class{DataFrame} of internal row metadata,
#' with number of rows equal to \code{nrow(x)}.
#' This is analogous to the user-visible \code{\link{colData}}.}
#' \item{\code{int_colData(x)}:}{Returns a \linkS4class{DataFrame} of internal column metadata,
#' with number of rows equal to \code{ncol(x)}.
#' This is analogous to the user-visible \code{\link{rowData}}.}
#' \item{\code{int_metadata(x)}:}{Returns a list of internal metadata, analogous to the user-visible \code{\link{metadata}}.}
#' } 
#'
#' @section Setters:
#' In the following code snippets, \code{x} is a \linkS4class{SingleCellExperiment}.
#' \describe{
#' \item{\code{int_elementMetadata(x) <- value}:}{Replaces the internal row metadata with \code{value},
#' a \linkS4class{DataFrame} with number of rows equal to \code{nrow(x)}.
#' This is analogous to the user-visible \code{\link{colData<-}}.}
#' \item{\code{int_colData(x) <- value}:}{Replaces the internal column metadata with \code{value},
#' a \linkS4class{DataFrame} with number of rows equal to \code{ncol(x)}.
#' This is analogous to the user-visible \code{\link{rowData<-}}.}
#' \item{\code{int_metadata(x) <- value}:}{Replaces the internal metadata with \code{value}, 
#' analogous to the user-visible \code{\link{metadata<-}}.}
#' }
#'
#' @section Comments:
#' The internal metadata fields allow easy and extensible storage of additional elements
#' that are parallel to the rows or columns of a \linkS4class{SingleCellExperiment} class.
#' This avoids the need to specify new slots and adjust the subsetting/combining code for a new data element.
#' For example, \code{\link{altExps}} and \code{\link{reducedDims}} are implemented as fields in the internal column metadata.
#'
#' That these elements are internal is important as this ensures that the implementation details are abstracted away.
#' Any user-interaction with these internal fields should be done via the designated getter and setter methods,
#' e.g., \code{\link{reducedDim}} and friends for retrieving or modifying reduced dimensions.
#' This provides developers with more freedom to change the internal representation without breaking user code.
#' 
#' Package developers intending to use these methods to store their own content should read the development vignette for guidance.
#' 
#' @seealso
#' \code{\link{colData}}, \code{\link{rowData}} and \code{\link{metadata}} for the user-visible equivalents.
#'
#' @author Aaron Lun
#'
#' @name internals
#' @docType methods
#' @aliases
#' int_colData
#' int_elementMetadata
#' int_metadata
#' int_colData,SingleCellExperiment-method
#' int_elementMetadata,SingleCellExperiment-method
#' int_metadata,SingleCellExperiment-method
#' int_colData<-
#' int_elementMetadata<-
#' int_metadata<-
#' int_colData<-,SingleCellExperiment-method
#' int_elementMetadata<-,SingleCellExperiment-method
#' int_metadata<-,SingleCellExperiment-method
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#' int_metadata(sce)$whee <- 1
NULL

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

