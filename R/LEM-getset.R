# Defines methods for new generics for LinearEmbeddingMatrix.

#' @export
setMethod("sampleFactors", "LinearEmbeddingMatrix", function(x) {
    x@sampleFactors
})

#' @export
setMethod("featureLoadings", "LinearEmbeddingMatrix", function(x, withDimnames=TRUE){
    fl <- x@featureLoadings
    if(withDimnames){
        colnames(fl) <- colnames(x@sampleFactors)
    }
    return(fl)
})

#' @export
setMethod("factorData", "LinearEmbeddingMatrix", function(x, withDimnames=TRUE){
    fd <- x@factorData
    if(withDimnames){
        rownames(fd) <- colnames(x@sampleFactors)
    }
    return(fd)
})

#' @export
setReplaceMethod("sampleFactors", "LinearEmbeddingMatrix", function(x, value) {
    x@sampleFactors <- value
    validObject(x)
    return(x)
})

#' @export
setReplaceMethod("featureLoadings", "LinearEmbeddingMatrix", function(x, value) {
    x@featureLoadings <- value
    validObject(x)
    return(x)
})

#' @export
setReplaceMethod("factorData", "LinearEmbeddingMatrix", function(x, value) {
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
    dimnames(sampleFactors(x))
})

#' @export
setReplaceMethod("dimnames", "LinearEmbeddingMatrix", function(x, value) {
    dimnames(sampleFactors(x)) <- value
    return(x)
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
