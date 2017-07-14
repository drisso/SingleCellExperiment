# Official getter/setter functions.

setMethod("reducedDims", "SingleCellExperiment", function(x) {
    x@reducedDims
})

setReplaceMethod("int_reducedDims", "SingleCellExperiment", function(x, value) {
    x@reducedDims <- value
    return(x)
})

setReplaceMethod("reducedDims", "SingleCellExperiment", function(x, value) {
    for (i in seq_along(value)) {
        if (!.not_reddim_mat(value[[i]], x)) { rownames(value[[i]]) <- colnames(x) }
    }
    int_reducedDims(x) <- value
    validObject(x)
    return(x)
})

setMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type) {
    reducedDims(x)[[type]]
})

setMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type) {
  reducedDims(x)[[type]]
})

setMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type) {
  r <- reducedDims(x)
  if(length(r) > 0) {
    return(r[[1]])
  } else {
    return(r)
  }
})

setReplaceMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    if (!.not_reddim_mat(value, x)) { rownames(value) <- colnames(x) }
    rd <- reducedDims(x)
    rd[[type]] <- value
    int_reducedDims(x) <- rd
    validObject(x)
    return(x)
})

# Internal getter/setter functions.

setMethod("int_elementMetadata", "SingleCellExperiment", function(x) x@int_elementMetadata)
setReplaceMethod("int_elementMetadata", "SingleCellExperiment", function(x, value) {
    x@int_elementMetadata <- value
    return(x)
})

setMethod("int_colData", "SingleCellExperiment", function(x) x@int_colData)
setReplaceMethod("int_colData", "SingleCellExperiment", function(x, value) {
    x@int_colData <- value
    return(x)
})

setMethod("int_metadata", "SingleCellExperiment", function(x) x@int_metadata)
setReplaceMethod("int_metadata", "SingleCellExperiment", function(x, value) {
    x@int_metadata <- value
    return(x)
})

# Size factor getter/setter functions.

setMethod("sizeFactors", "SingleCellExperiment", function(object, type=NULL) {
    field <- .get_sf_field(type)
    return(int_colData(object)[[field]])
})

setReplaceMethod("sizeFactors", "SingleCellExperiment", function(object, type=NULL, ..., value) {
    field <- .get_sf_field(type)
    cd <- int_colData(object)
    cd[[field]] <- value
    int_colData(object) <- cd
    return(object)
})

# Spike-in getter/setter functions.

setMethod("isSpike", c("SingleCellExperiment", "character"), function(x, type) {
    field <- .get_spike_field(type)
    if (!type %in% spikeNames(x)) {
        stop("spike-in set '", type, "' does not exist")
    }
    return(int_elementMetadata(x)[[field]])
})

for (sig in c("missing", "NULL")){
    setMethod("isSpike", c("SingleCellExperiment", sig), function(x, type) {
        return(int_elementMetadata(x)[[.spike_field]])
    })
}

setReplaceMethod("isSpike", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    md <- int_metadata(x)
    rd <- int_elementMetadata(x)
    field <- .get_spike_field(type)

    if (is.null(value)) {
        md$spike_names <- setdiff(md$spike_names, type) # Deleting if NULL.
        rd[[field]] <- NULL
    } else {
        md$spike_names <- union(md$spike_names, type)
        rd[[field]] <- .convert_subset_spike(value, .length=nrow(x), .names=rownames(x))
    }

    int_metadata(x) <- md
    int_elementMetadata(x) <- rd # need this, otherwise 'isSpike' below won't be up-to-date.

    # Updating is_spike.
    all.spikes <- lapply(spikeNames(x), isSpike, x=x)
    combined.spikes <- Reduce("|", all.spikes)
    rd[[.spike_field]] <- combined.spikes
    int_elementMetadata(x) <- rd
    return(x)
})

# colData / rowData

setMethod("colData", "SingleCellExperiment", function(x, internal=FALSE) {
  if(internal) {
    if (any(colnames(x@colData) %in% colnames(int_colData(x)))) {
      cn <- colnames(x@colData)[which(colnames(x@colData) %in%
                                        colnames(int_colData(x)))]
      warning("Overlapping column names (", paste(cn, collapse = ", "),
              ") between internal and external colData.")
    }
    cbind(callNextMethod(), int_colData(x))
  } else {
    callNextMethod()
  }
})

setMethod("rowData", "SingleCellExperiment", function(x, internal=FALSE) {
  if(internal) {
    if (any(colnames(mcols(x)) %in% colnames(int_elementMetadata(x)))) {
      cn <- colnames(mcols(x))[which(colnames(mcols(x)) %in%
                                                colnames(int_elementMetadata(x)))]
      warning("Overlapping column names (", paste(cn, collapse = ", "),
              ") between internal and external rowData.")
    }
    cbind(callNextMethod(), int_elementMetadata(x))
  } else {
    callNextMethod()
  }
})

# Other useful functions.

setMethod("spikeNames", "SingleCellExperiment", function(x) {
    int_metadata(x)$spike_names
})

setMethod("objectVersion", "SingleCellExperiment", function(x) {
    int_metadata(x)$version
})

setAs("SummarizedExperiment",
      "SingleCellExperiment",
      function(from) {
        SingleCellExperiment(assays = assays(from),
                             colData = colData(from),
                             rowData = rowData(from),
                             metadata = metadata(from)
                             )
      }
)
