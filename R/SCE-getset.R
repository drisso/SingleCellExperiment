#############################################
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

#############################################
# Other internal getter/setter functions.

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

#############################################
# Size factor getter/setter functions.

#' @export
#' @importFrom BiocGenerics sizeFactors
setMethod("sizeFactors", "SingleCellExperiment", function(object, type=NULL) {
    field <- .get_sf_field(type)
    return(int_colData(object)[[field]])
})

#' @export
#' @importFrom BiocGenerics "sizeFactors<-"
setReplaceMethod("sizeFactors", "SingleCellExperiment", function(object, type=NULL, ..., value) {
    field <- .get_sf_field(type)
    cd <- int_colData(object)
    cd[[field]] <- value
    int_colData(object) <- cd

    if (!is.null(type)) { 
        md <- int_metadata(object)
        if (is.null(value)) { 
            md$size_factor_names <- setdiff(md$size_factor_names, type)
        } else {
            md$size_factor_names <- union(md$size_factor_names, type)
        }
        int_metadata(object) <- md
    }
    return(object)
})

#' @export
setMethod("clearSizeFactors", "SingleCellExperiment", function(object) {
    sizeFactors(object) <- NULL

    cd <- int_colData(object)
    for (sf in sizeFactorNames(object)) {
        field <- .get_sf_field(sf)
        cd[[field]] <- NULL        
    }
    int_colData(object) <- cd

    md <- int_metadata(object) 
    md$size_factor_names <- character(0)
    int_metadata(object) <- md
    return(object)
})

# Spike-in getter/setter functions.

#' @export
setMethod("isSpike", c("SingleCellExperiment", "character"), function(x, type) {
    field <- .get_spike_field(type)
    return(int_elementMetadata(x)[[field]])
})

for (sig in c("missing", "NULL")){
    setMethod("isSpike", c("SingleCellExperiment", sig), function(x, type) {
        return(int_elementMetadata(x)[[.spike_field]])
    })
}

#' @export
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

#' @export
setMethod("clearSpikes", "SingleCellExperiment", function(x) {
    spike.sets <- spikeNames(x)
    rd <- int_elementMetadata(x)
    for (s in spike.sets) {
        field <- .get_spike_field(s)
        rd[[field]] <- NULL
    }
    rd[[.spike_field]] <- NULL # Emptying out the reference as well.
    int_elementMetadata(x) <- rd
    
    md <- int_metadata(x)
    md$spike_names <- character(0)
    int_metadata(x) <- md
    return(x)
})

#############################################
# colData / rowData getters, with options for accessing internal fields.

#' @export
#' @importFrom SummarizedExperiment colData
setMethod("colData", "SingleCellExperiment", function(x, internal=FALSE) {
    if(internal) {
        cn <- colnames(x@colData) # need explicit slot reference to avoid recursive colData() calling.
        conflict <- cn %in% colnames(int_colData(x))
        if (any(conflict)) {
            cn <- cn[conflict]
            if (length(cn) > 2) {
                cn <- c(cn[seq_len(2)], "...")
            }
            warning("overlapping names in internal and external colData (", paste(cn, collapse = ", "), ")")
        }
        cbind(callNextMethod(), int_colData(x))
    } else {
        callNextMethod()
  }
})

#' @export
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment rowData
setMethod("rowData", "SingleCellExperiment", function(x, internal=FALSE) {
    if(internal) {
        cn <- colnames(mcols(x))
        conflict <- cn %in% colnames(int_elementMetadata(x))
        if (any(conflict)) { 
            cn <- cn[conflict]
            if (length(cn) > 2) {
                cn <- c(cn[seq_len(2)], "...")
            }
            warning("overlapping names in internal and external rowData (", paste(cn, collapse = ", "), ")")
        }
        rv <- cbind(callNextMethod(), int_elementMetadata(x))
    } else {
        rv <- callNextMethod()
    }
    rownames(rv) <- rownames(x)
    rv
})

#############################################
# Other useful functions.

#' @export
setMethod("spikeNames", "SingleCellExperiment", function(x) {
    int_metadata(x)$spike_names
})

#' @export
setMethod("sizeFactorNames", "SingleCellExperiment", function(object) {
    int_metadata(object)$size_factor_names
})

#' @export
setMethod("objectVersion", "SingleCellExperiment", function(x) {
    int_metadata(x)$version
})

#############################################
# This defines some convenience wrappers for common entires in the assays slot.

GET_FUN <- function(exprs_values) {
    (exprs_values) # To ensure evaluation
    function(object) {
        assay(object, i=exprs_values)
    }
}

SET_FUN <- function(exprs_values) {
    (exprs_values) # To ensure evaluation
    function(object, value) {
        assay(object, i=exprs_values) <- value
        object
    }
}

#' @export
#' @importFrom BiocGenerics counts
setMethod("counts", "SingleCellExperiment", GET_FUN("counts"))

#' @export
#' @importFrom BiocGenerics "counts<-"
setReplaceMethod("counts", c("SingleCellExperiment", "ANY"), SET_FUN("counts"))

#' @export
setMethod("logcounts", "SingleCellExperiment", GET_FUN("logcounts"))

#' @export
setReplaceMethod("logcounts", c("SingleCellExperiment", "ANY"), SET_FUN("logcounts"))

#' @export
setMethod("normcounts", "SingleCellExperiment", GET_FUN("normcounts"))

#' @export
setReplaceMethod("normcounts", c("SingleCellExperiment", "ANY"), SET_FUN("normcounts"))

#' @export
setMethod("cpm", "SingleCellExperiment", GET_FUN("cpm"))

#' @export
setReplaceMethod("cpm", c("SingleCellExperiment", "ANY"), SET_FUN("cpm"))

#' @export
setMethod("tpm", "SingleCellExperiment", GET_FUN("tpm"))

#' @export
setReplaceMethod("tpm", c("SingleCellExperiment", "ANY"), SET_FUN("tpm"))
