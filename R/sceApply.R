#' Applying over parts of a SingleCellExperiment
#'
#' Apply a function over different pieces of data extracted from a \linkS4class{SingleCellExperiment}.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param which A list of \linkS4class{SCEInput} objects.
#' Each one may have different arguments to pass to \code{FUN}.
#' @param FUN A function to apply to each piece of data extracted with \code{\link{getInput}}.
#' @param ... Further (named) arguments to pass to \code{FUN}.
#' 
#' @return A list of the same length and names as \code{which},
#' containing the output of \code{FUN} applied to each corresponding piece of data.
#' 
#' @author Aaron Lun
#'
#' @details
#' If any named arguments in \code{...} overlap with arguments in a given entry \code{which}, the latter are used.
#' This enables \code{FUN} to use data-specific arguments while retaining the convenience of a single argument for globally applicable parameters.
#'
#' The \linkS4class{SCEInput} subclasses that can be used in \code{which} depend on what is supported by \code{FUN}.
#' For example, if a function can accept both matrix and SingleCellExperiment inputs, it can be used with, e.g., both \linkS4class{MainExpInput} and \linkS4class{AssayInput} objects.
#'
#' @examples
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' pca <- matrix(runif(ncells*5), ncells)
#' sce <- SingleCellExperiment(assays=list(counts=u),
#'     reducedDims=SimpleList(PCA=pca))
#' sce.copy <- SingleCellExperiment(assays=list(counts=u*10),
#'     reducedDims=SimpleList(PCA=pca/10))
#' altExp(sce, "BLAH") <- sce.copy
#'
#' # Here, using a very simple function that just 
#' # computes the mean of the input for each cell.
#' FUN <- function(y, transposed=FALSE, multiplier=1) {
#'     if (transposed) rowMeans(y) * multiplier
#'     else colMeans(y) * multiplier
#' }
#' 
#' # Applying over all of the specified parts of 'sce'.
#' sceApply(sce, FUN=FUN, which=list(
#'     AssayInput(assay=1),
#'     ReducedDimInput(type=1, transposed=TRUE),
#'     AltAssayInput(experiment=1, assay="counts", multiplier=2),
#'     AltReducedDimInput(experiment="BLAH", type="PCA", transposed=TRUE)
#' ), multiplier=1)
#' 
#' @seealso
#' \linkS4class{SCEInput} and the related class hierarchy, to specify different parts of the SingleCellExperiment.
#' @export
sceApply <- function(x, which, FUN, ...) {
    output <- which
    common <- list(...)

    for (i in seq_along(which)) {
        current <- which[[i]]
        args <- c(inputArguments(current), common)
        args <- args[!duplicated(names(args))]
        all.args <- c(list(getInput(x, current)), args)
        output[[i]] <- do.call(FUN, all.args)
    }

    output
}
