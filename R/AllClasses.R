#' @export
#' @importFrom utils packageVersion
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("SingleCellExperiment",
    slots=c(int_elementMetadata = "DataFrame",
        int_colData = "DataFrame",
        int_metadata = "list",
        reducedDims = "SimpleList"),
    contains = "RangedSummarizedExperiment",
    prototype = prototype(int_metadata=list(version=packageVersion("SingleCellExperiment"),
        spike_names=character(0),
        size_factor_names=character(0),
        reducedDims=setNames(SimpleList(), character(0)))
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
