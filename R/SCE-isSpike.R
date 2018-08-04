# Spike-in getter/setter functions.

.spike_field <- "is_spike"

.get_spike_field <- function(type, check=TRUE) {
    if (check && length(type)!=1L) {
        stop("'type' must be a character vector of length 1")
    }
    sprintf("%s_%s", .spike_field, type)
}

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
        output <- logical(nrow(x))
        names(output) <- rownames(x)
        output[value] <- TRUE
        rd[[field]] <- unname(output)
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

#' @export
setMethod("spikeNames", "SingleCellExperiment", function(x) {
    int_metadata(x)$spike_names
})
