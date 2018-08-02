.sce_validity <- function(object) {
    msg <- NULL

    # Checking dimensions of reduced coordinates.
    rd <- reducedDims(object, withDimnames=FALSE)
    if (length(rd)) {
        if (any(unlist(lapply(rd, .not_reddim_mat, object=object)))) {
            msg <- c(msg, "each element of 'reducedDims' must be a matrix-like object with nrow equal to 'ncol(object)'")
        }
    }

    # Checking dimensions of internal objects.
    if (nrow(int_elementMetadata(object))!=nrow(object)) {
        msg <- c(msg, "'nrow' of internal 'rowData' not equal to 'nrow(object)'")
    }
    if (nrow(int_colData(object))!=ncol(object)) {
        msg <- c(msg, "'nrow' of internal 'colData' not equal to 'ncol(object)'")
    }

    # Checking spike-in names are present and accounted for.
    spike.fields <- .get_spike_field(spikeNames(object), check=FALSE)
    lost.spikes <- ! spike.fields %in% colnames(int_elementMetadata(object))
    if (any(lost.spikes)) {
        was.lost <- spikeNames(object)[lost.spikes][1]
        msg <- c(msg, sprintf("no field specifying rows belonging to spike-in set '%s'", was.lost))
    }

    # Checking the size factor names as well.
    sf.fields <- vapply(sizeFactorNames(object), .get_sf_field, FUN.VALUE="")
    lost.sfs <- ! sf.fields %in% colnames(int_colData(object))
    if (any(lost.sfs)) {
        was.lost <- sizeFactorNames(object)[lost.sfs][1]
        msg <- c(msg, sprintf("no field specifying size factors for set '%s'", was.lost))
    }

    # Checking version.
    if (objectVersion(object) < "0.98.0") {
        msg <- c(msg, "object is out of date, update with 'updateSCE(object)'")
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("SingleCellExperiment", .sce_validity)
