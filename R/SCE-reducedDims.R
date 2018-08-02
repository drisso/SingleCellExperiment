# Getter/setter functions for reducedDims.

#' @export
setMethod("reducedDims", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    value <- x@reducedDims
    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- colnames(x)
        }
    }
    return(value)
})

setReplaceMethod("int_reducedDims", "SingleCellExperiment", function(x, value) {
    x@reducedDims <- value
    return(x)
})

.not_reddim_mat <- function(val, object) {
    return(is.null(val) || NROW(val)!=NCOL(object));
}

#' @export
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
setReplaceMethod("reducedDims", "SingleCellExperiment", function(x, value) {
    value <- as(value, "List")
    if (is.null(names(value))) {
        names(value) <- character(length(value))
    }

    int_reducedDims(x) <- value
    validObject(x)
    return(x)
})

#' @export
setMethod("reducedDimNames", "SingleCellExperiment", function(x) {
    names(reducedDims(x))
})

#' @export
setReplaceMethod("reducedDimNames", c("SingleCellExperiment", "character"), function(x, value) {
    out <- reducedDims(x, withDimnames=FALSE)
    names(out) <- value
    int_reducedDims(x) <- out
    return(x)
})

#' @export
setMethod("reducedDim", "SingleCellExperiment", function(x, type=1, withDimnames=TRUE) {
    r <- reducedDims(x, withDimnames=FALSE)
    if (length(r)==0) { 
        return(NULL)
    }

    out <- r[[type]]
    if (withDimnames && !is.null(out)) {
        rownames(out) <- colnames(x)
    }
    return(out)
})

#' @export
setReplaceMethod("reducedDim", "SingleCellExperiment", function(x, type=1, ..., value) {
    rd <- reducedDims(x, withDimnames=FALSE)
    if (is.numeric(type) && type > length(rd)+1) { 
        stop("subscript is out of bounds") 
    }
    rd[[type]] <- value
    int_reducedDims(x) <- rd
    validObject(x)
    return(x)
})
