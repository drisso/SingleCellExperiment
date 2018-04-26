#' @export
#' @importFrom utils packageVersion
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
                                                 size_factor_names=character(0)))
)

#' @export
#' @importClassesFrom S4Vectors DataFrame Annotated 
setClass("LinearEmbeddingMatrix",
         slots = c(sampleFactors = "ANY",
                   featureLoadings = "ANY",
                   factorData = "DataFrame"),
         contains = "Annotated")
