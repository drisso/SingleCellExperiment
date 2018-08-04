# colData / rowData getters, with options for accessing internal fields.

#' @export
#' @importFrom SummarizedExperiment colData
setMethod("colData", "SingleCellExperiment", function(x, ..., internal=FALSE) {
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
        cbind(callNextMethod(x, ...), int_colData(x))
    } else {
        callNextMethod(x, ...)
  }
})

#' @export
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment rowData
setMethod("rowData", "SingleCellExperiment", function(x, ..., internal=FALSE) {
    if (internal) {
        cn <- colnames(mcols(x))
        conflict <- cn %in% colnames(int_elementMetadata(x))
        if (any(conflict)) {
            cn <- cn[conflict]
            if (length(cn) > 2) {
                cn <- c(cn[seq_len(2)], "...")
            }
            warning("overlapping names in internal and external rowData (", paste(cn, collapse = ", "), ")")
        }
        cbind(callNextMethod(x, ...), int_elementMetadata(x))
    } else {
        callNextMethod(x, ...)
    }
})
