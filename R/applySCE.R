#' Applying over parts of a SingleCellExperiment
#'
#' Apply a function over different pieces of data extracted from a \linkS4class{SingleCellExperiment}.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param which A list of \linkS4class{SCEInput} objects.
#' Each entry may contain different arguments to pass to \code{FUN}.
#' @param FUN A function to apply to each piece of data extracted with \code{\link{getInput}}.
#' @param ... Further (named) arguments to pass to \code{FUN}.
#' @param SIMPLIFY Logical scalar indicating whether the output should be simplified.
#' 
#' @return 
#' In most cases or when \code{SIMPLIFY=FALSE}, a list of the same length and names as \code{which} is returned,
#' containing the output of \code{FUN} applied to each corresponding piece of data.
#' 
#' If \code{SIMPLIFY=TRUE} and certain conditions are fulfilled,
#' a SingleCellExperiment is returned where the results for each \code{which} are mapped to the relevant main or alternative Experiments of the output.
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
#' The default of \code{SIMPLIFY=TRUE} is intended as a user-level convenience when \code{FUN} returns a SingleCellExperiment with the same number of columns for all \code{which},
#' and \code{which} itself contains no more than one reference to the main or each alternative Experiment in \code{x}.
#' When these conditions are fulfilled, the results are collated into a single SingleCellExperiment for easier downstream manipulation.
#'
#' @section Developer note:
#' When using this function inside other functions, developers should set \code{SIMPLIFY=FALSE} to guarantee consistent output for arbitrary \code{which}.
#' If simplification is necessary, the output of this function can be explicitly passed to \code{\link{simplifyToSCE}},
#' typically with \code{\link{warn.level=3}} to throw an error if simplification is not possible.
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
#' applySCE(sce, FUN=FUN, which=list(
#'     AssayInput(assay=1),
#'     ReducedDimInput(type=1, transposed=TRUE),
#'     AltAssayInput(experiment=1, assay="counts", multiplier=2),
#'     AltReducedDimInput(experiment="BLAH", type="PCA", transposed=TRUE)
#' ), multiplier=1)
#' 
#' @seealso
#' \linkS4class{SCEInput} and the related class hierarchy, to specify different parts of the SingleCellExperiment.
#'
#' \code{\link{simplifyToSCE}}, which is used when \code{SIMPLIFY=TRUE}.
#' @export
applySCE <- function(x, which, FUN, ..., SIMPLIFY=TRUE) {
    output <- which
    common <- list(...)

    for (i in seq_along(which)) {
        current <- which[[i]]
        args <- c(inputArguments(current), common)
        args <- args[!duplicated(names(args))]
        all.args <- c(list(getInput(x, current)), args)
        output[[i]] <- do.call(FUN, all.args)
    }

    if (SIMPLIFY) {
        attempt <- simplifyToSCE(output, which, x, warn.level=1)
        if (!is.null(attempt)) {
            output <- attempt
        }
    }

    output
}
