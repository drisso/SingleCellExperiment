#' @export
#' @importFrom S4Vectors SimpleList 
#' @importFrom methods is as
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
SingleCellExperiment <- function(..., reducedDims=SimpleList()) {
    se <- SummarizedExperiment(...)
    if(!is(se, "RangedSummarizedExperiment")) {
        se <- as(se, "RangedSummarizedExperiment")
    }
    .rse_to_sce(se, reducedDims)
}

#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
#' @importFrom BiocGenerics nrow ncol
.rse_to_sce <- function(rse, reducedDims) {
    out <- new("SingleCellExperiment", rse, 
        int_elementMetadata=new("DataFrame", nrows=nrow(rse)),
        int_colData=new("DataFrame", nrows=ncol(rse)))
    reducedDims(out) <- reducedDims
    out
}

#' @exportMethod coerce
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom S4Vectors SimpleList
setAs("RangedSummarizedExperiment", "SingleCellExperiment", function(from) {
    .rse_to_sce(from, SimpleList())
})

#' @exportMethod coerce
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
setAs("SummarizedExperiment", "SingleCellExperiment", function(from) {
    .rse_to_sce(as(from, "RangedSummarizedExperiment"), SimpleList())
})
