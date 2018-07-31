# Getter/setter functions for reducedDims.

#' @export
setMethod("reducedDims", "SingleCellExperiment", function(x) {
    x@reducedDims
})

setReplaceMethod("int_reducedDims", "SingleCellExperiment", function(x, value) {
    x@reducedDims <- value
    return(x)
})

#' @export
setReplaceMethod("reducedDims", "SingleCellExperiment", function(x, value) {
    for (i in seq_along(value)) {
        if (!.not_reddim_mat(value[[i]], x)) { rownames(value[[i]]) <- colnames(x) }
    }
    int_reducedDims(x) <- value
    validObject(x)
    return(x)
})

#' @export
setMethod("reducedDimNames", "SingleCellExperiment", function(x) {
    rdn <- names(reducedDims(x))
    if (is.null(rdn)) {
        rdn <- character(length(reducedDims(x)))
    }
    return(rdn)
})

#' @export
setMethod("reducedDim", "SingleCellExperiment", function(x, type=1) {
    r <- reducedDims(x)
    if(length(r)==0) { return(NULL) }
    return(r[[type]])
})

#' @export
setReplaceMethod("reducedDim", "SingleCellExperiment", function(x, type=1, ..., value) {
    if (!.not_reddim_mat(value, x)) { rownames(value) <- colnames(x) }
    rd <- reducedDims(x)
    if (is.numeric(type) && type > length(rd)+1) { stop("subscript is out of bounds") }
    rd[[type]] <- value
    int_reducedDims(x) <- rd
    validObject(x)
    return(x)
})
