.alt_key <- "altExps"

#############################
# Getters.

#' @export
setMethod("altExpNames", "SingleCellExperiment", function(x) {
    if (!.alt_key %in% colnames(int_colData(x))) {
        character(0)
    }  else {
        colnames(int_colData(x)[[.alt_key]]) 
    }
})

#' @export
#' @importFrom S4Vectors List
#' @importClassesFrom S4Vectors SimpleList
setMethod("altExps", "SingleCellExperiment", function(x, withColData=TRUE) {
    if (!.alt_key %in% colnames(int_colData(x))) {
        List()
    }  else {
        out <- lapply(int_colData(x)[[.alt_key]], .get_se)
        if (withColData) {
            for (i in seq_along(out)) {
                colData(out[[i]]) <- colData(x)
            }
        }
        as(out, "SimpleList")
    }
})

#' @export
setMethod("altExp", "SingleCellExperiment", function(x, e=1, withColData=TRUE) {
    internals <- int_colData(x)
    if (!.alt_key %in% colnames(internals)) {
        stop("no alternative experiments in 'x'")
    }
    out <- .get_se(internals[,.alt_key][,e])
    if (withColData) {
        colData(out) <- colData(x)
    }
    colnames(out) <- colnames(x)
    out
})

#############################
# Setters.

#' @export
setReplaceMethod("altExpNames", "SingleCellExperiment", function(x, value) {
    if (!.alt_key %in% colnames(int_colData(x))) {
        if (length(value) > 0L) {
            stop("no alternative experiments in 'x' to rename")
        }
    }  else {
        colnames(int_colData(x)[[.alt_key]]) <- as.character(value)
    }
    x
})

#' @export
#' @importClassesFrom S4Vectors SimpleList
setReplaceMethod("altExps", "SingleCellExperiment", function(x, value) {
    collected <- int_colData(x)[,0]
    for (i in seq_along(value)) {
        collected[[i]] <- SummarizedExperimentByColumn(value[[i]])
    }
    if (!is.null(names(value))) {
        colnames(collected) <- names(value)
    } else {
        colnames(collected) <- character(length(value))
    }
    int_colData(x)[[.alt_key]] <- collected
    x
})

#' @export
setReplaceMethod("altExp", "SingleCellExperiment", function(x, e=1, ..., value) {
    internals <- int_colData(x)
    if (!is.null(value)) {
        value <- SummarizedExperimentByColumn(value)
    }
    if (!.alt_key %in% colnames(internals)) {
        internals[[.alt_key]] <- internals[,0]
    }
    internals[[.alt_key]][[e]] <- value
    int_colData(x) <- internals
    x
})
