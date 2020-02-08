#' Size factor methods
#' 
#' Gets or sets the size factors for all cells in a \linkS4class{SingleCellExperiment} object.
#' 
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param value A numeric vector of length equal to \code{ncol(object)}, containing size factors for all cells.
#' @param ... Additional arguments, currently ignored.
#' @param onAbsence String indicating an additional action to take when size factors are absent:
#' nothing (\code{"none"}), a warning (\code{"warn"}) or an error (\code{"error"}).
#' 
#' @details
#' A size factor is a scaling factor used to divide the raw counts of a particular cell to obtain normalized expression values,
#' thus allowing downstream comparisons between cells that are not affected by differences in library size or total RNA content.
#' The \code{sizeFactors} methods can be used to get or set size factors for all cells in a SingleCellExperiment object.
#'
#' When setting size factors, the values are stored in the \code{\link{colData}} as the \code{\link{sizeFactors}} field.
#' This name is chosen for general consistency with other packages (e.g., \pkg{DESeq2})
#' and to allow the size factors to be easily extracted from the \code{\link{colData}} for use as covariates.
#'
#' For developers, \code{onAbsence} is provided to make it easier to mandate that \code{object} has size factors.
#' This avoids silent \code{NULL} values that flow to the rest of the function and make debugging difficult.
#' 
#' @return
#' For \code{sizeFactors}, a numeric vector is returned containing size factors for all cells.
#' If no size factors are available, a \code{NULL} is returned (and/or a warning or error, depending on \code{onAbsence}).
#' 
#' For \code{sizeFactors<-}, a modified \code{object} is returned with size factors in its \code{\link{colData}}.
#' 
#' @seealso
#' \linkS4class{SingleCellExperiment}, for the underlying class definition. 
#'
#' \code{librarySizeFactors} from the \pkg{scater} package
#' or \code{computeSumFactors} from the \pkg{scran} package,
#' as examples of functions that compute the size factors.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#' sizeFactors(sce) <- runif(ncol(sce))
#' sizeFactors(sce)
#'
#' @name sizeFactors
NULL

.sf_field <- "sizeFactor"

#' @export
#' @rdname sizeFactors
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocGenerics sizeFactors
setMethod("sizeFactors", "SingleCellExperiment", function(object, onAbsence="none") {
    object <- updateObject(object)
    output <- colData(object)[[.sf_field]]
    .absent_action(object, val=output, fun="sizeFactors", onAbsence=onAbsence) 
    output
})

#' @export
#' @rdname sizeFactors
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom BiocGenerics sizeFactors<-
setReplaceMethod("sizeFactors", "SingleCellExperiment", function(object, ..., value) {
    object <- updateObject(object)
    colData(object)[[.sf_field]] <- value
    object
})
