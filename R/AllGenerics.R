# Getter/setters for reducedDim.

setGeneric("reducedDim", function(x, ...) standardGeneric("reducedDim"))
setGeneric("reducedDim<-", function(x, ..., value) standardGeneric("reducedDim<-"))
setGeneric("reducedDimNames", function(x) standardGeneric("reducedDimNames"))

setGeneric("reducedDims", function(x) standardGeneric("reducedDims"))
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

setGeneric("isSpike", function(x, type, ...) standardGeneric("isSpike"))
setGeneric("isSpike<-", function(x, type, ..., value) standardGeneric("isSpike<-"))

setGeneric("objectVersion", function(x) standardGeneric("objectVersion"))
setGeneric("spikeNames", function(x) standardGeneric("spikeNames"))

# Convenience assay getter/setters.

setGeneric("normcounts", function(object, ...) standardGeneric("normcounts"))
setGeneric("normcounts<-", function(object, ..., value) standardGeneric("normcounts<-"))
setGeneric("logcounts", function(object, ...) standardGeneric("logcounts"))
setGeneric("logcounts<-", function(object, ..., value) standardGeneric("logcounts<-"))

setGeneric("cpm", function(object, ...) standardGeneric("cpm"))
setGeneric("cpm<-", function(object, ..., value) standardGeneric("cpm<-"))
setGeneric("tpm", function(object, ...) standardGeneric("tpm"))
setGeneric("tpm<-", function(object, ..., value) standardGeneric("tpm<-"))

