.alt_key <- "altExps"

#############################
# Getters.

#' @export
setMethod("altExpNames", "SingleCellExperiment", function(x) {
    colnames(int_colData(x)[[.alt_key]])
})

#' @export
#' @importFrom S4Vectors List
#' @importClassesFrom S4Vectors SimpleList
setMethod("altExps", "SingleCellExperiment", function(x, withColData=TRUE) {
    out <- lapply(int_colData(x)[[.alt_key]], .get_se)
    if (withColData) {
        for (i in seq_along(out)) {
            colData(out[[i]]) <- colData(x)
        }
    }
    as(out, "SimpleList")
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "missing"), function(x, e, withColData=TRUE) {
    altExp(x, 1, withColData)
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "numeric"), function(x, e=1, withColData=TRUE) {
    internals <- int_colData(x)
    out <- .get_se(internals[,.alt_key][,e])
    if (withColData) {
        colData(out) <- colData(x)
    }
    colnames(out) <- colnames(x)
    out
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "character"), function(x, e, withColData=TRUE) {
    internals <- int_colData(x)
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
    colnames(int_colData(x)[[.alt_key]]) <- as.character(value)
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
setReplaceMethod("altExp", c("SingleCellExperiment", "missing"), function(x, e, ..., value) {
    altExp(x, 1L) <- value
    x
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "numeric"), function(x, e=1, ..., value) {
    internals <- int_colData(x)
    if (e[1] > ncol(internals[[.alt_key]])) {
        stop("invalid subscript 'type'\nsubscript out of bounds")
    }
    if (!is.null(value)) {
        value <- SummarizedExperimentByColumn(value)
    }
    internals[[.alt_key]][[e]] <- value
    int_colData(x) <- internals
    x
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "character"), function(x, e, ..., value) {
    internals <- int_colData(x)
    if (!is.null(value)) {
        value <- SummarizedExperimentByColumn(value)
    }
    internals[[.alt_key]][[e]] <- value
    int_colData(x) <- internals
    x
})
