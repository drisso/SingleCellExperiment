# Defines methods for new generics for LinearEmbeddingMatrix.

#' @export
setMethod("sampleFactors", "LinearEmbeddingMatrix", function(x, withDimnames=TRUE) {
    sf <- x@sampleFactors
    if (withDimnames) {
        dimnames(sf) <- dimnames(x)
    }
    return(sf)
})

#' @export
setMethod("featureLoadings", "LinearEmbeddingMatrix", function(x, withDimnames=TRUE){
    fl <- x@featureLoadings
    if(withDimnames){
        colnames(fl) <- colnames(x)
    }
    return(fl)
})

#' @export
setMethod("factorData", "LinearEmbeddingMatrix", function(x){
    x@factorData
})

#' @export
setReplaceMethod("sampleFactors", "LinearEmbeddingMatrix", function(x, value) {
    if (length(dim(value))!=2L) {
        stop("'value' should be a matrix-like object") 
    }
    if ((!is.null(colnames(value)) && !identical(colnames(value), colnames(x))) 
        || (!is.null(rownames(value)) && !identical(rownames(value), rownames(x))) ) {
        stop("'dimnames(value)' must match 'dimnames(x)' when setting sampleFactors")
    }

    x@sampleFactors <- value
    validObject(x)
    return(x)
})

#' @export
setReplaceMethod("featureLoadings", "LinearEmbeddingMatrix", function(x, value) {
    if (length(dim(value))!=2L) {
        stop("'value' should be a matrix-like object")
    }
    if (!is.null(colnames(value)) && !identical(colnames(value), colnames(x))){
        stop("'colnames(value)' must match 'colnames(x)' when setting featureLoadings")
    }
  
    x@featureLoadings <- value
    validObject(x)
    return(x)
})

#' @export
setReplaceMethod("factorData", "LinearEmbeddingMatrix", function(x, value) {
    if (!is(value, "DataFrame")) {
        stop("'value' should be a DataFrame")
    }
    if (!is.null(rownames(value)) && !identical(rownames(value), colnames(x))){
        stop("'rownames(value)' must match 'colnames(x)' when setting factorData")
    }
    
    x@factorData <- value
    validObject(x)
    return(x)
})

#' @export
setMethod("$", "LinearEmbeddingMatrix", function(x, name) {
    factorData(x)[[name]]
})

#' @export
setReplaceMethod("$", "LinearEmbeddingMatrix", function(x, name, value) {
    factorData(x)[[name]] <- value
    return(x)
})

#############################################
# Define matrix methods.

#' @export
setMethod("dim", "LinearEmbeddingMatrix", function(x) {
    dim(sampleFactors(x))
})

#' @export
setMethod("dimnames", "LinearEmbeddingMatrix", function(x) {
    list(x@NAMES, rownames(factorData(x)))
})

#' @export
setReplaceMethod("dimnames", "LinearEmbeddingMatrix", function(x, value) {
    fd <- factorData(x)
    rownames(fd) <- value[[2]]
    BiocGenerics:::replaceSlots(x, NAMES=value[[1]], factorData=fd, check=FALSE)
})

#' @export
#' @method as.matrix LinearEmbeddingMatrix
as.matrix.LinearEmbeddingMatrix <- function(x, ...) {
    y <- sampleFactors(x)
    while (!is.matrix(y)) { # in case of sparsity or other madness.
        y <- as.matrix(y)
    }
    return(y)
}

#' @export
setMethod("as.matrix", "LinearEmbeddingMatrix", as.matrix.LinearEmbeddingMatrix)
