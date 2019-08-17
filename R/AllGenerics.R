########################################
# Getters/setters for linearEmbedding.

#' @export
setGeneric("sampleFactors", function(x, ...) standardGeneric("sampleFactors"))

#' @export
setGeneric("sampleFactors<-", function(x, ..., value) standardGeneric("sampleFactors<-"))

#' @export
setGeneric("featureLoadings", function(x, ..., withDimnames=TRUE) standardGeneric("featureLoadings"))

#' @export
setGeneric("featureLoadings<-", function(x, ..., value) standardGeneric("featureLoadings<-"))

#' @export
setGeneric("factorData", function(x, ..., withDimnames=TRUE) standardGeneric("factorData"))

#' @export
setGeneric("factorData<-", function(x, ..., value) standardGeneric("factorData<-"))

########################################
# Getter/setters for reducedDim.

#' @export
setGeneric("reducedDim", function(x, i, ...) standardGeneric("reducedDim"))

#' @export
setGeneric("reducedDim<-", function(x, ..., value) standardGeneric("reducedDim<-"))

#' @export
setGeneric("reducedDimNames", function(x) standardGeneric("reducedDimNames"))

#' @export
setGeneric("reducedDimNames<-", function(x, value) standardGeneric("reducedDimNames<-"))

#' @export
setGeneric("reducedDims", function(x, ...) standardGeneric("reducedDims"))

#' @export
setGeneric("reducedDims<-", function(x, value) standardGeneric("reducedDims<-"))

setGeneric("int_reducedDims<-", function(x, value) standardGeneric("int_reducedDims<-"))

########################################
# Hidden getter/setters for internal slots.

#' @export
setGeneric("int_elementMetadata", function(x) standardGeneric("int_elementMetadata"))

#' @export
setGeneric("int_elementMetadata<-", function(x, value) standardGeneric("int_elementMetadata<-"))

#' @export
setGeneric("int_colData", function(x) standardGeneric("int_colData"))

#' @export
setGeneric("int_colData<-", function(x, value) standardGeneric("int_colData<-"))

#' @export
setGeneric("int_metadata", function(x) standardGeneric("int_metadata"))

#' @export
setGeneric("int_metadata<-", function(x, value) standardGeneric("int_metadata<-"))

########################################
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

########################################
# Alternative experiment getter/setters.

#' @export
setGeneric("altExp", function(x, ...) standardGeneric("altExp"))

#' @export
setGeneric("altExps", function(x, ...) standardGeneric("altExps"))

#' @export
setGeneric("altExpNames", function(x, ...) standardGeneric("altExpNames"))

#' @export
setGeneric("altExp<-", function(x, ..., value) standardGeneric("altExp<-"))

#' @export
setGeneric("altExps<-", function(x, ..., value) standardGeneric("altExps<-"))

#' @export
setGeneric("altExpNames<-", function(x, ..., value) standardGeneric("altExpNames<-"))

########################################
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

#' @export
setGeneric("weights<-", function(object, ..., value) standardGeneric("weights<-"))
