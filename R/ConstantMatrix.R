#' A DelayedMatrix containing a constant value
#'
#' A DelayedMatrix backend to efficiently mimic a matrix containing a constant value, without actually creating said matrix in memory.
#'
#' @param dim Integer vector of length 2, containing the dimensions of the matrix.
#' @param value Atomic scalar containing the value to fill the matrix.
#' @param seed A ConstantMatrixSeed object.
#' 
#' @return The \code{ConstantMatrixSeed} constructor returns a ConstantMatrixSeed object.
#'
#' The \code{ConstantMatrix} and \code{DelayedArray} constructors return a ConstantMatrix object.
#'
#' @details
#' This class allows us to efficiently create placeholder matrices, primarily to be used in \code{\link{unsplitAltExps}}.
#' For example, we can create a matrix full of \code{NA} values to enable combining of matrices without counterparts in other alternative experiments.
#'
#' We create this class instead of using the \linkS4class{RleArray},
#' as the latter requires some workarounds when the product of the dimensions is greater than the maximum integer value.
#'
#' @author Aaron Lun
#'
#' @examples
#' # This would ordinarily take up 8 TB of memory:
#' out <- ConstantMatrix(c(1e6, 1e6), value=NA_real_)
#' out
#'
#' @name ConstantMatrix
#' @aliases
#' ConstantMatrix-class
#' ConstantMatrixSeed-class
#' extract_array,ConstantMatrixSeed-method
NULL

#' @export
setClass("ConstantMatrixSeed", slots=c(dim='integer', value="ANY"))

setValidity("ConstantMatrixSeed", function(object) {
    msg <- character(0)

    if (length(dim(object)) != 2 || any(dim(object) < 0)) {
        msg <- c(msg, "'dim' must be an integer vector of length 2 with non-negative values")
    }

    if (length(object@value) != 1 || !is.atomic(object@value)) {
        msg <- c(msg, "'value' must be a single atomic value")
    }

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ConstantMatrix",
    contains="DelayedMatrix",
    representation(seed="ConstantMatrixSeed")
)

#' @export
#' @rdname ConstantMatrix
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "ConstantMatrixSeed", function(seed) new_DelayedArray(seed, Class="ConstantMatrix"))

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "ConstantMatrixSeed", function(x, index) {
    nr <- nrow(x)
    if (!is.null(index[[1]])) {
        nr <- length(index[[1]])
    }

    nc <- ncol(x)
    if (!is.null(index[[2]])) {
        nc <- length(index[[2]])
    }

    matrix(x@value, nr, nc) 
})

#' @export
#' @rdname ConstantMatrix
ConstantMatrixSeed <- function(dim, value) {
    new("ConstantMatrixSeed", dim=as.integer(dim), value=value)
}

#' @export
#' @rdname ConstantMatrix
ConstantMatrix <- function(dim, value) {
    DelayedArray(ConstantMatrixSeed(dim, value))
}
