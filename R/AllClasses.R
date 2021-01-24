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

#' @export
setClass("SCEInput", contains="VIRTUAL", slots=c(arguments="list"))

#' @export
setClass("MainExpInput", contains="SCEInput")

#' @export
setClass("AssayInput", contains="SCEInput", slots=c(assay="ANY"))

#' @export
setClass("ReducedDimInput", contains="SCEInput", slots=c(type="ANY"))

#' @export
setClass("AltExpInput", contains="SCEInput", slots=c(experiment="ANY"))

#' @export
setClass("AltAssayInput", contains="AltExpInput", slots=c(assay="ANY"))

#' @export
setClass("AltReducedDimInput", contains="AltExpInput", slots=c(type="ANY"))
