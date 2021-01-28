#' Applying over parts of a SingleCellExperiment
#'
#' Apply a function over the main and alternative Experiments of a \linkS4class{SingleCellExperiment}.
#'
#' @param X A \linkS4class{SingleCellExperiment} object.
#' @param WHICH A character or integer vector containing the names or positions of alternative Experiments to loop over.
#' @param FUN A function to apply to each Experiment.
#' @param ... Further (named) arguments to pass to all calls to \code{FUN}.
#' @param MAIN.ARGS A named list of arguments to pass to \code{FUN} for the main Experiment only.
#' Alternatively \code{NULL}, in which case the function is \emph{not} applied to the main Experiment.
#' @param ALT.ARGS A named list where each entry is named after an alternative Experiment and contains named arguments to use in \code{FUN} for that Experiment.
#' @param SIMPLIFY Logical scalar indicating whether the output should be simplified to a single SingleCellExperiment.
#' 
#' @return 
#' In most cases or when \code{SIMPLIFY=FALSE}, a list is returned containing the output of \code{FUN} applied to each Experiment.
#' If \code{MAIN.ARGS} is not \code{NULL}, the first entry corresponds to the result generated from the main Experiment;
#' all other results are generated according to the entries specified in \code{WHICH} and are named accordingly.
#' 
#' If \code{SIMPLIFY=TRUE} and certain conditions are fulfilled,
#' a SingleCellExperiment is returned where the results of \code{FUN} are mapped to the relevant main or alternative Experiments.
#' This mirrors the organization of Experiments in \code{X}.
#'
#' @author Aaron Lun
#'
#' @details
#' The behavior of this function is equivalent to creating a list containing \code{X} as the first entry and \code{\link{altExps}(X)} in the subsequent entries,
#' and then \code{\link{lapply}}ing over this list with \code{FUN} and the specified arguments.
#' In this manner, users can easily apply the same function to all the Experiments (main and alternative) in a \linkS4class{SingleCellExperiment} object.
#' 
#' Arguments in \code{...} are passed to all calls to \code{FUN}.
#' Arguments in \code{MAIN.ARGS} are only used in the call to \code{FUN} on the main Experiment.
#' Arguments in \code{ALT.ARGS} are passed to the call to \code{FUN} on the alternative Experiment of the same name.
#' For the last two, any arguments therein will override arguments of the same name in \code{...}.
#'
#' By default, looping is performed over all alternative Experiments, but the order and identities can be changed by setting \code{WHICH}.
#' Values of \code{WHICH} should be unique if any simplification of the output is desired.
#' If \code{MAIN.ARGS=NULL}, the main Experiment is ignored and the function is only applied to the alternative Experiments.
#'
#' The default of \code{SIMPLIFY=TRUE} is intended as a user-level convenience when all calls to \code{FUN} return a SingleCellExperiment with the same number of columns,
#' and \code{WHICH} itself contains no more than one reference to each alternative Experiment in \code{x}.
#' Under these conditions, the results are collated into a single SingleCellExperiment for easier downstream manipulation.
#'
#' @section Developer note:
#' When using this function inside other functions, developers should set \code{SIMPLIFY=FALSE} to guarantee consistent output for arbitrary \code{WHICH}.
#' If simplification is necessary, the output of this function can be explicitly passed to \code{\link{simplifyToSCE}},
#' typically with \code{warn.level=3} to throw an appropriate error if simplification is not possible.
#'
#' @examples
#' ncells <- 10
#' u <- matrix(rpois(200, 5), ncol=ncells)
#' sce <- SingleCellExperiment(assays=list(counts=u))
#' altExp(sce, "BLAH") <- SingleCellExperiment(assays=list(counts=u*10))
#' altExp(sce, "WHEE") <- SingleCellExperiment(assays=list(counts=u/10))
#'
#' # Here, using a very simple function that just 
#' # computes the mean of the input for each cell.
#' FUN <- function(y, multiplier=1) {
#'     colMeans(assay(y)) * multiplier
#' }
#' 
#' # Applying over all of the specified parts of 'sce'.
#' applySCE(sce, FUN=FUN)
#'
#' # Adding general arguments.
#' applySCE(sce, FUN=FUN, multiplier=5)
#' 
#' # Adding custom arguments.
#' applySCE(sce, FUN=FUN, MAIN.ARGS=list(multiplier=5))
#' applySCE(sce, FUN=FUN, ALT.ARGS=list(BLAH=list(multiplier=5)))
#'
#' # Skipping Experiments.
#' applySCE(sce, FUN=FUN, MAIN.ARGS=NULL) # skipping the main
#' applySCE(sce, FUN=FUN, WHICH=NULL) # skipping the alternatives
#' 
#' @seealso
#' \code{\link{simplifyToSCE}}, which is used when \code{SIMPLIFY=TRUE}.
#'
#' \code{\link{altExps}}, to manually extract the alternative Experiments for operations.
#' @export
applySCE <- function(X, FUN, WHICH=altExpNames(X), ..., MAIN.ARGS=list(), ALT.ARGS=list(), SIMPLIFY=TRUE) {
    use.main <- !is.null(MAIN.ARGS)
    output <- vector("list", use.main + length(WHICH))
    
    if (!is.character(WHICH)) {
        WHICH <- altExpNames(X)[WHICH]
    }
    ae.names <- WHICH
    if (use.main) {
        ae.names <- c("", WHICH)
    }
    names(output) <- ae.names

    common.args <- list(...)
    if (use.main) {
        tryCatch({
            output[[1]] <- do.call(FUN, c(list(X), .dedup_args(MAIN.ARGS, common.args)))
        }, error=function(err) {
            stop("'FUN' failed on the main Experiment:\n  ", conditionMessage(err))
        })
        ae.names <- c("", ae.names)
    }

    for (i in seq_along(WHICH)) {
        current <- WHICH[i]
        all.args <- c(list(altExp(X, current)), .dedup_args(ALT.ARGS[[current]], common.args))
        tryCatch({
            output[[i+use.main]] <- do.call(FUN, all.args)
        }, error=function(err) {
            stop("'FUN' failed on alternative Experiment '", current, "':\n  ", conditionMessage(err))
        })
    }

    if (SIMPLIFY) {
        which.main <- if (use.main) 1L else NULL
        attempt <- simplifyToSCE(output, which.main=which.main, warn.level=1)
        if (!is.null(attempt)) {
            output <- attempt
        }
    }

    output
}

.dedup_args <- function(arg1, arg2) {
    args <- c(arg1, arg2)
    args[!duplicated(names(args))]
}
