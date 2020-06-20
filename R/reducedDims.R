#' Reduced dimensions methods
#'
#' Methods to get or set dimensionality reduction results in a \linkS4class{SingleCellExperiment} object.
#' These are typically used to store and retrieve low-dimensional representations of single-cell datasets.
#' Each row of a reduced dimension result is expected to correspond to a column of the SingleCellExperiment object.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{reducedDim(x, type, withDimnames=TRUE)}:}{
#' Retrieves a matrix (or matrix-like object) containing reduced dimension coordinates for cells (rows) and dimensions (columns).
#' \code{type} is either a string specifying the name of the dimensionality reduction result in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired result.
#' If \code{withDimnames=TRUE}, row names of the output matrix are replaced with the column names of \code{x}.
#' }
#' \item{\code{reducedDimNames(x)}:}{
#' Returns a character vector containing the names of all dimensionality reduction results in \code{x}.
#' This is guaranteed to be of the same length as the number of results, though the names may not be unique.
#' }
#' \item{\code{reducedDims(x, withDimnames=TRUE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more dimensionality reduction results.
#' Each result is a matrix (or matrix-like object) with the same number of rows.
#' If \code{withDimnames=TRUE}, row names of each matrix are replaced with the column names of \code{x}.
#' }
#' }
#'
#' @section Single-result setter:
#' \code{reducedDim(x, type) <- value} will add or replace a dimensionality reduction result
#' in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{type} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{type} is missing, \code{value} is assigned to the first result.
#' If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed1"}.
#' \item If \code{type} is a numeric scalar, it must be within the range of existing results, and \code{value} will be assigned to the result at that index.
#' \item If \code{type} is a string and a result exists with this name, \code{value} is assigned to to that result.
#' Otherwise a new result with this name is append to the existing list of results.
#' }
#'
#' \code{value} is expected to be a matrix or matrix-like object with number of rows equal to \code{ncol(x)}.
#' Alternatively, if \code{value} is \code{NULL}, the result corresponding to \code{type} is removed from the object.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{reducedDims(x) <- value}:}{
#' Replaces all dimensionality reduction results in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of matrices or matrix-like objects
#' with number of rows equal to \code{ncol(x)}.
#'
#' If \code{value} is named, those names will be used to name the dimensionality reduction results in \code{x}.
#' Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed"}.
#'
#' If \code{value} is \code{NULL}, all dimensionality reduction results in \code{x} are removed.
#' }
#' \item{\code{reducedDimNames(x) <- value}:}{
#' Replaces all names for dimensionality reduction results in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of results currently in \code{x}.
#' }
#' }
#'
#' @author Aaron Lun and Kevin Rue-Albrecht
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE)
#' reducedDim(sce, "PCA")
#' reducedDim(sce, "tSNE")
#' reducedDims(sce)
#'
#' reducedDim(sce, "PCA") <- NULL
#' reducedDims(sce)
#'
#' reducedDims(sce) <- SimpleList()
#' reducedDims(sce)
#'
#' @name reducedDims
#' @docType methods
#' @aliases
#' reducedDim reducedDims reducedDimNames
#' reducedDim,SingleCellExperiment,missing-method
#' reducedDim,SingleCellExperiment,numeric-method
#' reducedDim,SingleCellExperiment,character-method
#' reducedDims,SingleCellExperiment-method
#' reducedDimNames,SingleCellExperiment-method
#' reducedDim<- reducedDims<- reducedDimNames<-
#' reducedDim<-,SingleCellExperiment,missing-method
#' reducedDim<-,SingleCellExperiment,numeric-method
#' reducedDim<-,SingleCellExperiment,character-method
#' reducedDims<-,SingleCellExperiment-method
#' reducedDimNames<-,SingleCellExperiment,character-method
NULL

.red_key <- "reducedDims"

#' @export
setMethod("reducedDims", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    value <- .get_internal_all(x, 
        getfun=int_colData, 
        key=.red_key)

    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- colnames(x)
        }
    }
    value
})

#' @export
setReplaceMethod("reducedDims", "SingleCellExperiment", function(x, value) {
    .set_internal_all(x, value, 
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.red_key,
        convertfun=NULL,
        xdimfun=ncol,
        vdimfun=nrow,
        funstr="reducedDims",
        xdimstr="ncol",
        vdimstr="rows")
})

#' @export
setMethod("reducedDimNames", "SingleCellExperiment", function(x) {
    .get_internal_names(x, 
        getfun=int_colData, 
        key=.red_key)
})

#' @export
setReplaceMethod("reducedDimNames", c("SingleCellExperiment", "character"), function(x, value) {
    .set_internal_names(x, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.red_key)
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE) {
    .get_internal_missing(x, 
        basefun=reducedDim, 
        namefun=reducedDimNames, 
        funstr="reducedDim",
        withDimnames=withDimnames)
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE) {
    out <- .get_internal_integer(x, type,
        getfun=int_colData,
        key=.red_key,
        funstr="reducedDim",
        substr="type")

    if (withDimnames) {
        rownames(out) <- colnames(x)
    }

    out
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    out <- .get_internal_character(x, type,
        getfun=int_colData,
        key=.red_key,
        funstr="reducedDim",
        substr="type",
        namestr="reducedDimNames")

    if (withDimnames) {
        rownames(out) <- colnames(x)
    }

    out
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    .set_internal_missing(x, value,
        basefun=`reducedDim<-`,
        namefun=reducedDimNames
    )
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type, ..., value) {
    .set_internal_numeric(x, type, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.red_key,
        convertfun=NULL,
        xdimfun=ncol,
        vdimfun=nrow,
        funstr="reducedDim",
        xdimstr="ncol",
        vdimstr="rows",
        substr="type")
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    .set_internal_character(x, type, value, 
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.red_key,
        convertfun=NULL,
        xdimfun=ncol, 
        vdimfun=nrow,
        funstr="reducedDim", 
        xdimstr="ncol",
        vdimstr="rows", 
        substr="type") 
})
