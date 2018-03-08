# Getter/setters for reducedDim.

#' @export
setGeneric("reducedDim", function(x, ...) standardGeneric("reducedDim"))

#' @export
setGeneric("reducedDim<-", function(x, ..., value) standardGeneric("reducedDim<-"))

#' @export
setGeneric("reducedDimNames", function(x) standardGeneric("reducedDimNames"))

#' @export
setGeneric("reducedDims", function(x) standardGeneric("reducedDims"))

#' @export
setGeneric("reducedDims<-", function(x, value) standardGeneric("reducedDims<-"))

setGeneric("int_reducedDims<-", function(x, value) standardGeneric("int_reducedDims<-"))

# Hidden getter/setters for internal slots.

setGeneric("int_elementMetadata", function(x) standardGeneric("int_elementMetadata"))

setGeneric("int_elementMetadata<-", function(x, value) standardGeneric("int_elementMetadata<-"))

setGeneric("int_colData", function(x) standardGeneric("int_colData"))

setGeneric("int_colData<-", function(x, value) standardGeneric("int_colData<-"))

setGeneric("int_metadata", function(x) standardGeneric("int_metadata"))

setGeneric("int_metadata<-", function(x, value) standardGeneric("int_metadata<-"))

# Loose getter/setters (i.e., with no official slot).

#' @export
setGeneric("isSpike", function(x, type, ...) standardGeneric("isSpike"))

#' @export
setGeneric("isSpike<-", function(x, type, ..., value) standardGeneric("isSpike<-"))

#' @export
setGeneric("clearSpikes", function(x, type, ..., value) standardGeneric("clearSpikes"))

#' @export
setGeneric("clearSizeFactors", function(object) standardGeneric("clearSizeFactors"))

#' @export
setGeneric("objectVersion", function(x) standardGeneric("objectVersion"))

#' @export
setGeneric("spikeNames", function(x) standardGeneric("spikeNames"))

#' @export
setGeneric("sizeFactorNames", function(object) standardGeneric("sizeFactorNames"))

# Convenience assay getter/setters.

#' @export
setGeneric("normcounts", function(object, ...) standardGeneric("normcounts"))

#' @export
setGeneric("normcounts<-", function(object, ..., value) standardGeneric("normcounts<-"))

#' @export
setGeneric("logcounts", function(object, ...) standardGeneric("logcounts"))

#' @export
setGeneric("logcounts<-", function(object, ..., value) standardGeneric("logcounts<-"))

#' @export
setGeneric("cpm", function(object, ...) standardGeneric("cpm"))

#' @export
setGeneric("cpm<-", function(object, ..., value) standardGeneric("cpm<-"))

#' @export
setGeneric("tpm", function(object, ...) standardGeneric("tpm"))

#' @export
setGeneric("tpm<-", function(object, ..., value) standardGeneric("tpm<-"))
