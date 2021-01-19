#' Swap main and alternative Experiments
#'
#' Swap the main Experiment for an alternative Experiment in a \linkS4class{SingleCellExperiment} object.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param name String or integer scalar specifying the alternative Experiment to use to replace the main Experiment.
#' @param saved String specifying the name to use to save the original \code{x} as an alternative experiment in the output.
#' If \code{NULL}, the original is not saved.
#' @param withColData Logical scalar specifying whether the column metadata of \code{x} should be preserved in the output.
#'
#' @return A SingleCellExperiment derived from \code{altExp(x, name)}.
#' This contains all alternative Experiments in \code{altExps(x)}, excluding the one that was promoted to the main Experiment.
#' An additional alternative Experiment containing \code{x} may be included if \code{saved} is specified.
#' 
#' @details
#' During the course of an analysis, we may need to perform operations on each of the alternative Experiments in turn.
#' This would require us to repeatedly call \code{altExp(x, name)} prior to running downstream functions on those Experiments.
#' In such cases, it may be more convenient to switch the main Experiment with the desired alternative Experiments,
#' allowing a particular section of the analysis to be performed on the latter by default.
#' 
#' For example, the initial phases of the analysis might use the entire set of features.
#' At some point, we might want to focus only on a subset of features of interest,
#' but we do not want to discard the rest of the features.
#' This can be achieved by storing the subset as an alternative Experiment and swapping it with the main Experiment,
#' as shown in the Examples below.
#'
#' If \code{withColData=TRUE}, the column metadata of the output object is set to \code{colData(x)}.
#' As a side-effect, any column data previously \code{altExp(x, name)} is stored in the \code{saved} alternative Experiment of the output.
#' This is necessary to preserve the column metadata while achieving reversibility (see below).
#' Setting \code{withColData=FALSE} will omit the \code{colData} exchange.
#'
#' \code{swapAltExp} is almost perfectly reversible, i.e., \code{swapAltExp(swapAltExp(x, name, saved), saved, name)} should return something very similar to \code{x}. 
#' The only exceptions are that the order of \code{\link{altExpNames}} is changed,
#' and that any non-\code{NULL} \code{\link{mainExpName}} in \code{altExp(x, name)} will be lost.
#'
#' @author Aaron Lun
#' @seealso 
#' \code{\link{altExps}}, for a description of the alternative Experiment concept.
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # using the class example
#'
#' # Let's say we defined a subset of genes of interest.
#' # We can save the feature set as its own altExp.
#' hvgs <- 1:10
#' altExp(sce, "subset") <- sce[hvgs,] 
#'
#' # At some point, we want to do our analysis on the HVGs only,
#' # but we want to hold onto the other features for later reference.
#' sce <- swapAltExp(sce, name="subset", saved="all")
#' sce
#' 
#' # Once we're done, it is straightforward to switch back.
#' swapAltExp(sce, "all") 
#' @export
swapAltExp <- function(x, name, saved=mainExpName(x), withColData=TRUE) {
    y <- altExp(x, name, withColData=withColData)
    y <- as(y, "SingleCellExperiment")

    all.ae <- altExps(x, withColData=FALSE)
    altExps(y) <- all.ae[setdiff(names(all.ae), name)]

    if (!is.null(saved)) {
        old <- removeAltExps(x)
        if (withColData) {
            colData(old) <- colData(altExp(x, name))
        }
        mainExpName(old) <- NULL
        altExp(y, saved) <- old
    }

    mainExpName(y) <- name
    y
}
