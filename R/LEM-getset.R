# Defines methods for new generics for LinearEmbeddingMatrix.

#' @export
setMethod("sampleFactors", "LinearEmbeddingMatrix", function(x) {
    x@sampleFactors
})

#' @export
setMethod("featureLoadings", "LinearEmbeddingMatrix", function(x){
    x@featureLoadings
})

#' @export
setMethod("factorData", "LinearEmbeddingMatrix", function(x){
    x@factorData
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
