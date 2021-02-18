#' The reduced.dim.matrix class
#'
#' A matrix class that retains its attributes upon being subsetted or combined.
#' This is useful for storing metadata about a dimensionality reduction result alongside the matrix,
#' and for ensuring that the metadata persists when the matrix is stored inside \code{\link{reducedDims}}.
#'
#' @section Constructor:
#' \code{reduced.dim.matrix(x, ...)} will return a reduced.dim.matrix object, given a matrix input \code{x}.
#' Arguments in \code{...} should be named and are stored as custom attributes in the output.
#' Any arguments named \code{dim} or \code{dimnames} are ignored.
#'
#' @section Subsetting:
#' \code{x[i, j, ..., drop=FALSE]} will subset a reduced.dim.matrix \code{x} in the same manner as a base matrix.
#' The only difference is that a reduced.dim.matrix will be returned, retaining any custom attributes in \code{x}.
#' Note that no custom attributes are retained if the return value is a vector with \code{drop=TRUE}.
#'
#' @section Combining:
#' \code{rbind(...)} will combine multiple reduced.dim.matrix inputs in \code{...} by row, 
#' while \code{cbind(...)} will combine those inputs by column.
#' 
#' If the custom attributes are the same across all objects \code{...},
#' a reduced.dim.matrix is returned containing all combined rows/columns as well as the custom attributes.
#'
#' If the custom attributes are different, a warning is issued.
#' A matrix is returned containing all combined rows/columns; no custom attributes are retained.
#' 
#' @author Aaron Lun
#'
#' @examples
#' # Typical PC result, with metadata stored in the attributes:
#' pc <- matrix(runif(500), ncol=5)
#' attr(pc, "sdev") <- 1:100
#' attr(pc, "rotation") <- matrix(rnorm(20), ncol=5)
#'
#' # Disappears upon subsetting and combining!
#' attributes(pc[1:10,])
#' attributes(rbind(pc, pc))
#'
#' # Transformed into a reduced.dim.matrix:
#' rd.pc <- reduced.dim.matrix(pc)
#'
#' attributes(rd.pc[1:10,])
#' attributes(rbind(rd.pc, rd.pc))
#'
#' @seealso
#' \code{\link{reducedDims}}, to store these objects in a \linkS4class{SingleCellExperiment}.
#'
#' @aliases
#' reduced.dim.matrix
#' reduced.dim.matrix-class
#' [.reduced.dim.matrix
#' cbind.reduced.dim.matrix
#' rbind.reduced.dim.matrix
#'
#' @name reduced.dim.matrix
NULL

#' @export
reduced.dim.matrix <- function(x, ...) {
    class(x) <- c("reduced.dim.matrix", "matrix")
    mostattributes(x) <- c(attributes(x), list(...))
    x
}

#' @rawNamespace exportClasses(reduced.dim.matrix)
setOldClass(c("reduced.dim.matrix", "matrix"))

#' @export
`[.reduced.dim.matrix` <- function(x, i, j, ..., drop=FALSE) {
    at <- attributes(x)
    out <- NextMethod()
    if (!is.null(dim(out))) {
        at <- at[setdiff(names(at), c("dim", "dimnames"))]
        mostattributes(out) <- c(attributes(out), at)
    }
    out
}

.check_reddim_attributes <- function(available) {
    all.attr <- lapply(available, attributes)

    # Ignore dim and dimnames.
    for (i in seq_along(all.attr)) {
        current <- all.attr[[i]] 
        all.attr[[i]] <- current[setdiff(names(current), c("dim", "dimnames"))]
    }

    u.attr <- unique(all.attr)
    if (length(u.attr) > 1) {
        warning("mismatched custom attributes when combining 'reduced.dim.matrix' objects")
        NULL
    } else {
        u.attr[[1]]
    }
}

#' @export
#' @method rbind reduced.dim.matrix
rbind.reduced.dim.matrix <- function(..., deparse.level=1) {
    available <- list(...)
    u.attr <- .check_reddim_attributes(available)

    available <- lapply(available, unclass)
    out <- do.call(rbind, available)

    mostattributes(out) <- c(attributes(out), u.attr)
    out
}

#' @export
#' @method cbind reduced.dim.matrix
cbind.reduced.dim.matrix <- function(..., deparse.level=1) {
    available <- list(...)
    u.attr <- .check_reddim_attributes(available)

    available <- lapply(available, unclass)
    out <- do.call(cbind, available)

    mostattributes(out) <- c(attributes(out), u.attr)
    out
}
