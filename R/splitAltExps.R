#' Split off alternative features
#'
#' Split a \linkS4class{SingleCellExperiment} based on the feature type, creating alternative Experiments to hold features that are not in the majority set.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param f A character vector or factor of length equal to \code{nrow(x)}, specifying the feature type of each row.
#' @param ref String indicating which level of \code{f} should be treated as the main set.
#'
#' @return
#' A SingleCellExperiment where each row corresponds to a feature in the main set.
#' Each other feature type is stored as an alternative Experiment, accessible by \code{\link{altExp}}.
#' \code{ref} is used as the \code{\link{mainExpName}}.
#'
#' @details
#' This function provides a convenient way to create a SingleCellExperiment with alternative Experiments.
#' For example, a SingleCellExperiment with rows corresponding to all features can be quickly split into endogenous genes (main)
#' and other alternative features like spike-in transcripts and antibody tags.
#'
#' By default, the most frequent level of \code{f} is treated as the \code{ref} if the latter is not specified.
#'
#' @author
#' Aaron Lun
#'
#' @seealso
#' \code{\link{altExp}}, to access and manipulate the alternative Experiment fields.
#'
#' \code{\link{unsplitAltExps}}, to reverse the splitting.
#' @examples
#' example(SingleCellExperiment, echo=FALSE)
#' feat.type <- sample(c("endog", "ERCC", "CITE"), nrow(sce),
#'     replace=TRUE, p=c(0.8, 0.1, 0.1))
#'
#' sce2 <- splitAltExps(sce, feat.type)
#' sce2
#' @export
#' @importFrom SummarizedExperiment colData colData<-
splitAltExps <- function(x, f, ref=NULL) {
    by.feat <- split(seq_along(f), f)
    if (is.null(ref)) {
        ref <- names(by.feat)[which.max(lengths(by.feat))]
    }

    x0 <- x[by.feat[[ref]], ]
    for (other in setdiff(names(by.feat), ref)) {
        # Clearing out the colData() before adding it.
        subset <- x[by.feat[[other]],]
        colData(subset) <- colData(subset)[, 0]
        mainExpName(subset) <- NULL
        altExp(x0, other) <- subset
    }

    mainExpName(x0) <- ref
    x0
}
