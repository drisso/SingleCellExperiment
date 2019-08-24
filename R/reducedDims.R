#' Reduced dimensions methods
#'
#' Methods to get or set dimensionality reduction results in a \linkS4class{SingleCellExperiment} object.
#' 
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{reducedDim(x, type, withDimnames=TRUE)}:}{
#' Retrieves a numeric matrix containing reduced dimension coordinates for cells (rows) and dimensions (columns).
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
#' Each result is a numeric matrix with the same number of rows.
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
#' \code{value} is expected to be a numeric matrix-like object with number of rows equal to \code{ncol(x)}.
#' Alternatively, if \code{value} is \code{NULL}, the set of results corresponding to \code{type} is removed from the object.
#' 
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{reducedDims(x) <- value}:}{
#' Replaces all dimensionality reduction results in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of numeric matrix-like objects 
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
#' reducedDimNames<-,SingleCellExperiment-method
NULL

# Getter/setter functions for reducedDims.

.red_key <- "reducedDims"
.unnamed <- "unnamed"

#' @export
#' @importFrom S4Vectors List
#' @importClassesFrom S4Vectors SimpleList
setMethod("reducedDims", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    x <- updateObject(x)
    value <- as(int_colData(x)[[.red_key]], "SimpleList")
    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- colnames(x)
        }
    }
    value
})

#' @export
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame
setReplaceMethod("reducedDims", "SingleCellExperiment", function(x, value) {
    x <- updateObject(x)

    if (length(value)==0L) {
        collected <- int_colData(x)[,0]
    } else {
        nrows <- vapply(value, nrow, FUN.VALUE = 0L)
        if (!all(nrows == ncol(x))) {
            stop("elements of replacement 'reducedDims' do not have the correct number of rows")
        }
        collected <- do.call(DataFrame, lapply(value, I))
        if (is.null(names(value))) {
            colnames(collected) <- paste0(.unnamed, seq_along(value))
        }
    }

    int_colData(x)[[.red_key]] <- collected
    x
})

#' @export
setMethod("reducedDimNames", "SingleCellExperiment", function(x) {
    x <- updateObject(x)
    colnames(int_colData(x)[[.red_key]])
})

#' @export
setReplaceMethod("reducedDimNames", c("SingleCellExperiment", "character"), function(x, value) {
    x <- updateObject(x)
    colnames(int_colData(x)[[.red_key]]) <- value
    x
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE) {

    if (identical(length(reducedDimNames(x)), 0L)) {
        .Deprecated(msg="NULL is deprecated.")
        return(NULL)
        # To deprecate NULL and throw an error instead, remove the two lines above.
        stop(
            "'reducedDim(<", class(x), ">, ...) ",
            "length(reducedDims(<", class(x), ">)) is 0'")
    }

    reducedDim(x, 1L, withDimnames)
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type=1, withDimnames=TRUE) {
    x <- updateObject(x)
    internals <- int_colData(x)[[.red_key]]

    out <- tryCatch({
        internals[, type]
    }, error=function(err) {
        .Deprecated(msg="NULL is deprecated.")
        return(NULL)
        # To deprecate NULL and throw an error instead, remove the two lines above.
        stop("'reducedDim(<", class(x), ">, type=\"numeric\", ...)' ",
             "invalid subscript 'type'\n", conditionMessage(err))
    })

    if (withDimnames) {
        rownames(out) <- colnames(x)
    }

    out
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    msg <- paste0("'reducedDim(<", class(x), ">, type=\"character\", ...)' ",
                  "invalid subscript 'type'")

    x <- updateObject(x)
    internals <- int_colData(x)[[.red_key]]

    out <- tryCatch({
        internals[, type]
    }, error=function(err) {
        .Deprecated(msg="NULL is deprecated.")
        return(NULL)
        # To deprecate NULL and throw an error instead, remove the two lines above.
        stop(msg, "\n'", type, "' not in names(reducedDims(<", class(x), ">))")
    })

    if (withDimnames) {
        rownames(out) <- colnames(x)
    }

    out
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    type <- 1L

    if (0L == length(reducedDimNames(x))){
        ## Implementation 1
        # https://github.com/drisso/SingleCellExperiment/pull/35#issuecomment-522258649
        type <- paste0(.unnamed, 1L)
        ## Implementation 2
        # `reducedDims<-`` above creates character(n) names
        # type <- character(1)
        ## Implementation 3
        # This would implement the same behaviour as SummarizedExperiment::assay<-
        # stop("'reducedDim(<", class(x), ">) <- value' ", "length(reducedDims(<",
        #      class(x), ">)) is 0")
    }

    reducedDim(x, type) <- value
    x
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type = 1, ..., value) {
    x <- updateObject(x)

    if (!is.null(value) && !identical(nrow(value), ncol(x))) {
        stop("replacement 'reducedDim' has a different number of rows than 'ncol(x)'")
    }

    internals <- int_colData(x)

    if (type[1] > ncol(internals[[.red_key]])) {
        stop("invalid subscript 'type'\nsubscript out of bounds")
    }

    internals[[.red_key]][[type]] <- value
    int_colData(x) <- internals
    x
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    x <- updateObject(x)

    internals <- int_colData(x)
    if (!is.null(value) && !identical(nrow(value), ncol(x))) {
        stop("replacement 'reducedDim' has a different number of rows than 'ncol(x)'")
    }

    internals[[.red_key]][[type]] <- value
    int_colData(x) <- internals
    x
})
