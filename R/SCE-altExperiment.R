.alt_key <- "altExperiments"

.get_alt_experiment <- function(x, e) {
    internals <- int_colData(x)
    if (!.alt_key %in% colnames(internals)) {
        stop("no alternative experiments in 'x'")
    }
    out <- .get_se(internals[,.alt_key][,e])
    colnames(out) <- colnames(x)
    out
}

.set_alt_experiment <- function(x, e, se) {
    internals <- int_colData(x)
    y <- SummarizedExperimentByColumn(se)
    if (!.alt_key %in% colnames(internals)) {
        internals[[.alt_key]] <- internals[,0]
    }
    internals[[.alt_key]][[e]] <- y
    int_colData(x) <- internals
    x
}

#############################
# Getters.

#' @export
setMethod("altExperimentNames", "SingleCellExperiment", function(x) {
    if (!.alt_key %in% colnames(int_colData(x))) {
        character(0)
    }  else {
        colnames(int_colData(x)[[.alt_key]]) 
    }
})

#' @export
setMethod("altExperiment", "SingleCellExperiment", function(x, e=1) {
    .get_alt_experiment(x, e)
})

#' @export
#' @importFrom SummarizedExperiment assays
setMethod("altAssays", "SingleCellExperiment", function(x, e=1, ...) {
    assays(.get_alt_experiment(x, e), ...)
})

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("altAssay", "SingleCellExperiment", function(x, e=1, i=1, ...) {
    assay(.get_alt_experiment(x, e), i, ...)
})

#' @export
#' @importFrom SummarizedExperiment assayNames
setMethod("altAssayNames", "SingleCellExperiment", function(x, e=1, ...) {
    assayNames(.get_alt_experiment(x, e), ...)
})

#' @export
#' @importFrom SummarizedExperiment rowData
setMethod("altRowData", "SingleCellExperiment", function(x, e=1, ...) {
    rowData(.get_alt_experiment(x, e), ...)
})

#' @export
setMethod("altRowNames", "SingleCellExperiment", function(x, e=1) {
    rownames(.get_alt_experiment(x, e))
})

#############################
# Setters.

#' @export
setReplaceMethod("altExperimentNames", "SingleCellExperiment", function(x, value) {
    if (!.alt_key %in% colnames(int_colData(x)) && length(value) > 0L) {
        stop("no alternative experiments in 'x' to rename")
    } 
    colnames(int_colData(x)[[.alt_key]]) <- value
    x
})

#' @export
setReplaceMethod("altExperiment", "SingleCellExperiment", function(x, e=1, ..., value) 
    .set_alt_experiment(x, e, value)
)

#' @export
#' @importFrom SummarizedExperiment assays<-
setReplaceMethod("altAssays", "SingleCellExperiment", function(x, e=1, ..., value) {
    se <- .get_alt_experiment(x, e)
    assays(se, ...) <- value
    .set_alt_experiment(x, e, se)
})

#' @export
#' @importFrom SummarizedExperiment assay<-
setReplaceMethod("altAssay", "SingleCellExperiment", function(x, e=1, i=1, ..., value) {
    se <- .get_alt_experiment(x, e)    
    assay(se, i, ...) <- value
    .set_alt_experiment(x, e, se)
})

#' @export
#' @importFrom SummarizedExperiment assayNames<-
setReplaceMethod("altAssayNames", "SingleCellExperiment", function(x, e=1, ..., value) {
    se <- .get_alt_experiment(x, e)
    assayNames(se, ...) <- value
    .set_alt_experiment(x, e, se)
})

#' @export
#' @importFrom SummarizedExperiment rowData<-
setReplaceMethod("altRowData", "SingleCellExperiment", function(x, e=1, ..., value) {
    se <- .get_alt_experiment(x, e)    
    rowData(se, ...) <- value
    .set_alt_experiment(x, e, se)
})

#' @export
setReplaceMethod("altRowNames", "SingleCellExperiment", function(x, e=1, ..., value) {
    se <- .get_alt_experiment(x, e)    
    rownames(se) <- value
    .set_alt_experiment(x, e, se)
})
