#' @export
#' @importFrom S4Vectors SimpleList 
#' @importFrom methods is as
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
SingleCellExperiment <- function(..., reducedDims=SimpleList(), altExps=list()) {
    se <- SummarizedExperiment(...)
    if(!is(se, "RangedSummarizedExperiment")) {
        se <- as(se, "RangedSummarizedExperiment")
    }
    .rse_to_sce(se, reducedDims, altExps)
}

#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
#' @importFrom BiocGenerics nrow ncol
.rse_to_sce <- function(rse, reducedDims=SimpleList(), altExps=SimpleList()) {
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    
    out <- new("SingleCellExperiment", rse, 
        int_elementMetadata=new("DataFrame", nrows=nrow(rse)),
        int_colData=new("DataFrame", nrows=ncol(rse)))

    reducedDims(out) <- reducedDims
    altExps(out) <- altExps
    out
}

#' @exportMethod coerce
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setAs("RangedSummarizedExperiment", "SingleCellExperiment", function(from) {
    .rse_to_sce(from)
})

#' @exportMethod coerce
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SummarizedExperiment
setAs("SummarizedExperiment", "SingleCellExperiment", function(from) {
    .rse_to_sce(as(from, "RangedSummarizedExperiment"))
})
