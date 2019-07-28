#' @export
setMethod("[", c("SingleCellExperiment", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) {
        ii <- .convert_subset_index(i, rownames(x))
        int_elementMetadata(x) <- int_elementMetadata(x)[ii,,drop=FALSE]
    }

    if (!missing(j)) {
        jj <- .convert_subset_index(j, colnames(x))
        int_colData(x) <- int_colData(x)[jj,,drop=FALSE]
    }

    callNextMethod()
})

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData colData
setMethod("[<-", c("SingleCellExperiment", "ANY", "ANY", "SingleCellExperiment"), function(x, i, j, ..., value) {
    if (missing(i) && missing(j)) {
        return(value)
    }

    if (!missing(i)) {
        left <- .create_shell_rowdata(x)
        right <- .create_shell_rowdata(value)
        ii <- .convert_subset_index(i, rownames(x))

        tryCatch({
            left[ii,] <- right
        }, error=function(err) {
            stop("failed to replace 'int_elementMetadata'")
        })
        int_elementMetadata(x) <- rowData(left)
    }

    if (!missing(j)) {
        left <- .create_shell_coldata(x)
        right <- .create_shell_coldata(value)
        jj <- .convert_subset_index(j, colnames(x))

        tryCatch({
            left[,jj] <- right
        }, error=function(err) {
            stop("failed to replace 'int_colData'")
        })
        int_colData(x) <- colData(left)
    }

    int_metadata(x) <- int_metadata(value)
    callNextMethod()
})
