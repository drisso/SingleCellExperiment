#' @export
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
SingleCellExperiment <- function(..., reducedDims=SimpleList()) {
    se <- SummarizedExperiment(...)
    if(!is(se, "RangedSummarizedExperiment")) {
      rse <- as(se, "RangedSummarizedExperiment")
      rowData(rse) <- rowData(se)
    } else {
      rse <- se
    }
    out <- new("SingleCellExperiment", rse, 
        int_elementMetadata=DataFrame(matrix(0, nrow(se), 0)),
        int_colData=DataFrame(matrix(0, ncol(se), 0)))
    reducedDims(out) <- reducedDims
    return(out)
}

#' @exportMethod coerce
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors metadata
setAs("SummarizedExperiment", "SingleCellExperiment", function(from) {
    SingleCellExperiment(assays = assays(from),
                         colData = colData(from),
                         rowData = rowData(from),
                         metadata = metadata(from)
                         )
})
