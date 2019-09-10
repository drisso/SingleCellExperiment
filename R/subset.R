#' @export
setMethod("[", c("SingleCellExperiment", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    x <- updateObject(x)
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
    x <- updateObject(x)
    value <- updateObject(value)
    if (missing(i) && missing(j)) {
        return(value)
    }

    if (!missing(i)) {
        left <- int_elementMetadata(x)
        right <- int_elementMetadata(value)
        ii <- .convert_subset_index(i, rownames(x))

        tryCatch({
            left[ii,] <- right
        }, error=function(err) {
            stop(
                "failed to replace 'int_elementMetadata' in '<", class(x), ">[i,] <- value'\n",
                conditionMessage(err))
        })
        int_elementMetadata(x) <- left
    }

    if (!missing(j)) {
        left <- int_colData(x)
        right <- int_colData(value)
        jj <- .convert_subset_index(j, colnames(x))

        tryCatch({
            left[jj,] <- right
        }, error=function(err) {
            stop(
                "failed to replace 'int_colData' in '<", class(x), ">[,j] <- value'\n",
                conditionMessage(err))
        })
        int_colData(x) <- left
    }

    int_metadata(x) <- int_metadata(value)
    callNextMethod()
})
