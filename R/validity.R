.sce_validity <- function(object) {
    msg <- NULL

    if (nrow(int_elementMetadata(object))!=nrow(object)) {
        msg <- c(msg, "'nrow' of 'int_elementMetadata' not equal to 'nrow(object)'")
    }
    if (nrow(int_colData(object))!=ncol(object)) {
        msg <- c(msg, "'nrow' of 'int_colData' not equal to 'ncol(object)'")
    }

    if (objectVersion(object) >= "1.7.1") {
        if (!.red_key %in% colnames(int_colData(object))) {
            msg <- c(msg, "no 'reducedDims' field in 'int_colData'")
        }
        if (!.alt_key %in% colnames(int_colData(object))) {
            msg <- c(msg, "no 'altExps' field in 'int_colData'")
        }
    }

    if (objectVersion(object) >= "1.11.3") {
        if (!.rowp_key %in% colnames(int_elementMetadata(object))) {
            msg <- c(msg, "no 'rowPairs' field in 'int_elementMetadata'")
        }
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("SingleCellExperiment", .sce_validity)
