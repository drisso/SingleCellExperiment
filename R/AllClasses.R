#' @export
#' @rdname SingleCellExperiment
#' @importFrom utils packageVersion
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("SingleCellExperiment",
    slots=c(int_elementMetadata = "DataFrame",
        int_colData = "DataFrame",
        int_metadata = "list"),
    contains = "RangedSummarizedExperiment",
    prototype = prototype(
        int_metadata=list(
            version=packageVersion("SingleCellExperiment"),
            mainExpName=NULL
        )
    )
)

#' @export
#' @importClassesFrom S4Vectors DataFrame Annotated character_OR_NULL
setClass("LinearEmbeddingMatrix",
         slots = c(sampleFactors = "ANY",
                   featureLoadings = "ANY",
                   NAMES = "character_OR_NULL",
                   factorData = "DataFrame"),
         contains = "Annotated")

setClass("SummarizedExperimentByColumn", slots=c(se="SummarizedExperiment"))

setClass("DualSubset", slots=c(hits="SelfHits"))
