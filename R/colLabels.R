#' Get or set column labels
#'
#' Get or set column labels in an instance of a \linkS4class{SingleCellExperiment} class.
#' Labels are expected to represent information about the the biological state of each cell.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param value Any vector-like object of length equal to \code{ncol(object)}, containing labels for all cells.
#' Alternatively \code{NULL}, in which case existing label information is removed.
#' @param ... Additional arguments, currently ignored.
#' @param onAbsence String indicating an additional action to take when labels are absent:
#' nothing (\code{"none"}), a warning (\code{"warn"}) or an error (\code{"error"}).
#'
#' @details
#' A frequent task in single-cell data analyses is to label cells with some annotation,
#' e.g., cluster identities, predicted cell type classifications and so on.
#' In a \linkS4class{SummarizedExperiment}, the \code{\link{colData}} represents the ideal place for such annotations,
#' which can be easily set and retrieved with standard methods, e.g., \code{x$label <- my.labels}.
#'
#' That said, it is desirable to have some informal standardization of the name of the column used to store these annotations as this makes it easier to programmatically set sensible defaults for retrieval of the labels in downstream functions.
#' To this end, the \code{colLabels} function will get or set labels from the \code{"label"} field of the \code{\link{colData}}.
#' This considers the use case where there is a \dQuote{primary} set of labels that represents the default grouping of cells in downstream analyses.
#'
#' To illustrate, let's say we have a downstream function that accepts a SingleCellExperiment object and requires labels.
#' When defining our function, we can set \code{colLabels(x)} as the default value for our label argument.
#' This pattern is useful as it accommodates on-the-fly changes to a secondary set of labels in \code{x} without requiring the user to run \code{colLabels(x) <- second.labels}, while facilitating convenient use of the primary labels by default.
#' 
#' For developers, \code{onAbsence} is provided to make it easier to mandate that \code{x} actually has labels.
#' This avoids silent \code{NULL} values that flow to the rest of the function and make debugging difficult.
#'
#' @author Aaron Lun
#'
#' @return
#' For \code{colLabels}, a vector or equivalent is returned containing label assignments for all cells.
#' If no labels are available, a \code{NULL} is returned (and/or a warning or error, depending on \code{onAbsence}).
#' 
#' For \code{colLabels<-}, a modified \code{x} is returned with labels in its \code{\link{colData}}.
#' 
#' @seealso
#' \linkS4class{SingleCellExperiment}, for the underlying class definition. 
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#' colLabels(sce) <- sample(LETTERS, ncol(sce), replace=TRUE)
#' colLabels(sce)
#'
#' @docType methods
#' @name colLabels
#' @aliases colLabels colLabels<-
NULL

.label_field <- "label"

#' @export
#' @rdname colLabels
#' @importFrom SummarizedExperiment colData
setMethod("colLabels", "SingleCellExperiment", function(x, onAbsence="none") {
    output <- colData(x)[[.label_field]]
    .absent_action(x, val=output, fun="colLabels", onAbsence=onAbsence)
    output
})

#' @export
#' @rdname colLabels
#' @importFrom SummarizedExperiment colData<- colData
setReplaceMethod("colLabels", "SingleCellExperiment", function(x, ..., value) {
    colData(x)[[.label_field]] <- value
    x
})
