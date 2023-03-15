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
#' In the following code snippets, \code{x} is a SingleCellExperiment and \code{...} contains multiple SingleCellExperiment objects.
#' \describe{
#' \item{\code{combineCols(x, ..., delayed=TRUE, fill=NA, use.names=TRUE)}:}{
#' Returns a SingleCellExperiment where all objects are flexibly combined by column.
#' The assays and \code{\link{colData}} are combined as described in \code{?"\link{combineCols,SummarizedExperiment-method}"},
#' where assays or DataFrame columns missing in any given object are filled in with missing values before combining.
#'
#' Entries of the \code{\link{reducedDims}} with the same name across objects are combined by row.
#' If a dimensionality reduction result is not present for a particular SingleCellExperiment, it is represented by a matrix of \code{NA} values instead.
#' If corresponding \code{\link{reducedDim}} entries cannot be combined, e.g., due to inconsistent dimensions, they are omitted from the \code{\link{reducedDims}} of the output object with a warning.
#'
#' Entries of the \code{\link{altExps}} with the same name across objects are combined by column using the relevant \code{\link{combineCols}} method.
#' If a named entry is not present for a particular SingleCellExperiment, it is represented by a SummarizedExperiment with a single assay full of \code{fill} values.
#' If entries cannot be combined, e.g., due to inconsistent dimensions, they are omitted from the \code{\link{altExps}} of the output object with a warning.
#'
#' Entries of the \code{\link{colPairs}} with the same name across objects are concatenated together after adjusting the indices for each column's new position in the combined object.
#' If a named entry is not present for a particular SingleCellExperiments, it is assumed to contribute no column pairings and is ignored.
#' 
#' Entries of the \code{\link{rowPairs}} with the same name should be identical across objects if \code{use.names=FALSE}.
#' If \code{use.names=TRUE}, we attempt to merge together entries with the same name by taking the union of all column pairings.
#' However, if the same cell has a different set of pairings across objects, a warning is raised and we fall back to the \code{\link{rowPair}} entry from the first object.
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
#' combineCols,SingleCellExperiment-method
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
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    out <- callNextMethod()

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

    BiocGenerics:::replaceSlots(out, int_colData=int_cd, int_elementMetadata=int_em,
        int_metadata=int_m, check=FALSE)
})

#' @export
#' @importFrom BiocGenerics rbind cbind
setMethod("rbind", "SingleCellExperiment", function(..., deparse.level=1) {
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    out <- callNextMethod()

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

#' @export
#' @importFrom DelayedArray ConstantArray
#' @importFrom S4Vectors combineCols combineRows combineUniqueCols
setMethod("combineCols", "SingleCellExperiment", function(x, ..., delayed=TRUE, fill=NA, use.names=TRUE) {
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    ans <- callNextMethod()

    args <- list(x, ...)
    args <- lapply(args, updateObject)
    new.int_m <- do.call(c, unname(lapply(args, int_metadata)))

    # Merging everything but the alternative experiments and reduced dimensions.
    all.int_cd <- all.int_ed <- all.altexp <- all.reddim <- all.colp <- all.rowp <- vector("list", length(args))
    for (i in seq_along(all.int_cd)) {
        current.int_cd <- int_colData(args[[i]])
        current.int_ed <- int_elementMetadata(args[[i]])

        all.altexp[i] <- list(current.int_cd[[.alt_key]])
        all.reddim[i] <- list(current.int_cd[[.red_key]])
        all.colp[i] <- list(current.int_cd[[.colp_key]])
        all.rowp[i] <- list(current.int_ed[[.rowp_key]])

        current.int_cd[[.alt_key]] <- NULL
        current.int_cd[[.red_key]] <- NULL
        current.int_cd[[.colp_key]] <- NULL
        current.int_ed[[.rowp_key]] <- NULL

        all.int_cd[[i]] <- current.int_cd
        all.int_ed[[i]] <- current.int_ed
    }

    new.int_cd <- do.call(combineRows, all.int_cd)

    # Reduced dimensions need a bit more care.
    all.reddim.names <- Reduce(union, lapply(all.reddim, colnames))
    new.reddim <- new.int_cd[,0]

    for (i in all.reddim.names) {
        collated <- vector("list", length(all.reddim))
        first <- NULL

        for (j in seq_along(collated)) {
            if (i %in% names(all.reddim[[j]])) {
                collated[[j]] <- all.reddim[[j]][[i]]
                if (is.null(first)) {
                    first <- collated[[j]]
                }
            }
        }

        for (j in seq_along(collated)) {
            if (is.null(collated[[j]])) {
                collated[[j]] <- matrix(NA_real_, ncol(args[[j]]), ncol(first))
            }
        }

        tryCatch({
            new.reddim[[i]] <- do.call(rbind, collated)
        }, error=function(e) {
            warning("failed to combine '", i, "' in 'reducedDims(<", class(ans), ">)':\n  ", conditionMessage(e))
        })
    }

    new.int_cd[[.red_key]] <- new.reddim

    # Alternative experiments need a bit more care.
    all.altexp.names <- Reduce(union, lapply(all.altexp, colnames))
    new.altexps <- new.int_cd[,0]

    for (i in all.altexp.names) {
        collated <- vector("list", length(all.altexp))
        first <- NULL

        for (j in seq_along(collated)) {
            if (i %in% names(all.altexp[[j]])) {
                collated[[j]] <- all.altexp[[j]][[i]]@se
                if (is.null(first)) {
                    first <- collated[[j]]
                }
            }
        }

        for (j in seq_along(collated)) {
            if (is.null(collated[[j]])) {
                dummy.assay <- assays(first)[1]
                dummy.assay[[1]] <- ConstantArray(c(nrow(first), ncol(args[[j]])), fill)
                dummy <- SummarizedExperiment(dummy.assay, rowData=rowData(first)[,0], colData=colData(args[[j]])[,0])
                collated[[j]] <- as(dummy, class(first)) # probably not the best placeholder, but whatever.
            }
        }

        tryCatch({
            new.altexps[[i]] <- SummarizedExperimentByColumn(do.call(combineCols, c(collated, list(use.names=use.names, delayed=delayed, fill=fill))))
        }, error=function(e) {
            warning("failed to combine '", i, "' in 'altExps(<", class(ans), ">)':\n  ", conditionMessage(e))
        })
    }

    new.int_cd[[.alt_key]] <- new.altexps

    # As do the colPairs.
    all.colp.names <- Reduce(union, lapply(all.colp, colnames))
    new.colp <- new.int_cd[,0]

    for (i in all.colp.names) {
        collated <- vector("list", length(all.colp))
        first <- NULL

        for (j in seq_along(collated)) {
            if (i %in% names(all.colp[[j]])) {
                collated[[j]] <- all.colp[[j]][[i]]
            }
        }

        for (j in seq_along(collated)) {
            if (is.null(collated[[j]])) {
                collated[[j]] <- DualSubset(SelfHits(integer(0), integer(0), nnode=ncol(args[[j]])))
            }
        }

        tryCatch({
            new.colp[[i]] <- do.call(c, collated)
        }, error=function(e) {
            warning("failed to combine '", i, "' in 'colPairs(<", class(ans), ">)':\n  ", conditionMessage(e))
        })
    }

    new.int_cd[[.colp_key]] <- new.colp

    # As do the rowPairs.
    for (i in seq_along(all.rowp)) {
        rn <- rownames(args[[i]])
        rownames(all.int_ed[[i]]) <- rn
        rownames(all.rowp[[i]]) <- rn
    }
    new.rowp <- do.call(combineUniqueCols, c(all.rowp, list(use.names=use.names)))
    new.int_ed <- do.call(combineUniqueCols, c(all.int_ed, list(use.names=use.names)))
    rownames(new.rowp) <- rownames(new.int_ed) <- NULL
    new.int_ed[[.rowp_key]] <- new.rowp

    ans <- as(ans, "SingleCellExperiment")
    BiocGenerics:::replaceSlots(ans, int_colData=new.int_cd, int_elementMetadata=new.int_ed,
        int_metadata=new.int_m, check=FALSE)
})
