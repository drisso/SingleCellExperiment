#' @title Combining or subsetting SingleCellExperiment objects
#'
#' @description
#' An overview of methods to combine multiple \linkS4class{SingleCellExperiment} objects by row or column,
#' or to subset a SingleCellExperiment by row or column.
#' These methods are useful for ensuring that all data fields remain synchronized
#' when cells or genes are added or removed.
#'
#' @section Combining:
#' In the following code snippets, \code{...} contains one or more \linkS4class{SingleCellExperiment} objects.
#' \describe{
#' \item{\code{rbind(..., deparse.level=1)}:}{Returns a SingleCellExperiment where all objects in \code{...} are combined row-wise,
#' i.e., rows in successive objects are appended to the first object.
#'
#' Refer to \code{?"\link{rbind,SummarizedExperiment-method}"} for details on how metadata is combined in the output object.
#' Refer to \code{?\link[base]{rbind}} for the interpretation of \code{deparse.level}.
#'
#' Note that all objects in \code{...} must have the exact same values for \code{\link{reducedDims}} and \code{\link{altExps}}.
#' Any \code{\link{sizeFactors}} should either be \code{NULL} or contain the same values across objects.
#' }
#' \item{\code{cbind(..., deparse.level=1)}:}{Returns a SingleCellExperiment where
#' all objects in \code{...} are combined column-wise,
#' i.e., columns in successive objects are appended to the first object.
#'
#' Each object \code{x} in \code{...} must have the same values of \code{\link{reducedDimNames}(x)} (though they can be unordered).
#' Dimensionality reduction results with the same name across objects
#' will be combined row-wise to create the corresponding entry in the output object.
#'
#' Each object \code{x} in \code{...} must have the same values of \code{\link{altExpNames}(x)} (though they can be unordered).
#' Alternative Experiments with the same name across objects
#' will be combined column-wise to create the corresponding entry in the output object.
#'
#' \code{\link{sizeFactors}} should be either set to \code{NULL} in all objects, or set to a numeric vector in all objects.
#'
#' Refer to \code{?"\link{cbind,SummarizedExperiment-method}"} for details on how metadata is combined in the output object.
#' Refer to \code{?\link[base]{cbind}} for the interpretation of \code{deparse.level}.
#' }
#' }
#'
#' @section Subsetting:
#' In the following code snippets, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{x[i, j, ..., drop=TRUE]}:}{Returns a SingleCellExperiment containing the
#' specified rows \code{i} and columns \code{j}.
#'
#' \code{i} and \code{j} can be a logical, integer or character vector of subscripts,
#' indicating the rows and columns respectively to retain.
#' Either can be missing, in which case subsetting is only performed in the specified dimension.
#' If both are missing, no subsetting is performed.
#'
#' Arguments in \code{...} and \code{drop} are passed to to \code{\link{[,SummarizedExperiment-method}}.}
#' \item{\code{x[i, j, ...] <- value}:}{Replaces all data for rows \code{i} and columns {j}
#' with the corresponding fields in a SingleCellExperiment \code{value}.
#'
#' \code{i} and \code{j} can be a logical, integer or character vector of subscripts,
#' indicating the rows and columns respectively to replace.
#' Either can be missing, in which case replacement is only performed in the specified dimension.
#' If both are missing, \code{x} is replaced entirely with \code{value}.
#'
#' If \code{j} is specified, \code{value} is expected to have the same name and order of \code{\link{reducedDimNames}}
#' and \code{\link{altExpNames}} as \code{x}.
#' If \code{sizeFactors} is set for \code{x}, it should also be set for \code{value}.
#'
#' Arguments in \code{...} are passed to the corresponding \linkS4class{SummarizedExperiment} method.}
#' }
#'
#' @author
#' Aaron Lun
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # using the class example
#'
#' # Combining:
#' rbind(sce, sce)
#' cbind(sce, sce)
#'
#' # Subsetting:
#' sce[1:10,]
#' sce[,1:5]
#'
#' sce2 <- sce
#' sce2[1:10,] <- sce[11:20,]
#'
#' # Can also use subset()
#' sce$WHEE <- sample(LETTERS, ncol(sce), replace=TRUE)
#' subset(sce, , WHEE=="A")
#'
#' # Can also use split()
#' split(sce, sample(LETTERS, nrow(sce), replace=TRUE))
#'
#' @docType methods
#' @aliases
#' cbind,SingleCellExperiment-method
#' rbind,SingleCellExperiment-method
#' [,SingleCellExperiment,ANY-method
#' [,SingleCellExperiment,ANY,ANY-method
#' [,SingleCellExperiment,ANY,ANY,ANY-method
#' [<-,SingleCellExperiment,ANY,ANY,SingleCellExperiment-method
#'
#' @name SCE-combine
#' @rdname combine
NULL

#' @export
#' @importFrom BiocGenerics rbind cbind
setMethod("cbind", "SingleCellExperiment", function(..., deparse.level=1) {
    args <- list(...)
    args <- lapply(args, updateObject)
    int_m <- do.call(c, unname(lapply(args, int_metadata)))

    tryCatch({
        int_cd <- do.call(rbind, lapply(args, int_colData))
    }, error=function(err) {
        stop(
            "failed to combine 'int_colData' in 'cbind(<", class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })

    # Creating a shell to avoid having to pull out .cbind.DataFrame
    # to fuse metadata along the dimension not being combined.
    row_shells <- lapply(args, .create_shell_rowdata)
    tryCatch({
        combined <- do.call(cbind, row_shells)
    }, error=function(err) {
        stop(
            "failed to combine 'int_elementMetadata' in 'cbind(<", class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    int_em <- rowData(combined)

    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, int_colData=int_cd, int_elementMetadata=int_em,
        int_metadata=int_m, check=FALSE)
})

#' @export
#' @importFrom BiocGenerics rbind cbind
setMethod("rbind", "SingleCellExperiment", function(..., deparse.level=1) {
    args <- list(...)
    args <- lapply(args, updateObject)
    int_m <- do.call(c, unname(lapply(args, int_metadata)))

    tryCatch({
        int_em <- do.call(rbind, lapply(args, int_elementMetadata))
    }, error=function(err) {
        stop(
            "failed to combine 'int_elementMetadata' in 'rbind(<", class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })

    # Creating a shell to avoid having to pull out .cbind.DataFrame
    # to fuse metadata along the dimension not being combined.
    col_shells <- lapply(args, .create_shell_coldata)
    tryCatch({
        combined <- do.call(rbind, col_shells)
    }, error=function(err) {
        stop(
            "failed to combine 'int_colData' in 'rbind(<", class(args[[1]]), ">)'\n",
            conditionMessage(err))
    })
    int_cd <- colData(combined)

    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, int_colData=int_cd, int_elementMetadata=int_em,
        int_metadata=int_m, check=FALSE)
})

#' @importFrom SummarizedExperiment SummarizedExperiment
.create_shell_coldata <- function(x) {
    SummarizedExperiment(colData=int_colData(x))
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.create_shell_rowdata <- function(x) {
    SummarizedExperiment(rowData=int_elementMetadata(x))
}

