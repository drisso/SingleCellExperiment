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
setMethod("[<-", c("SingleCellExperiment", "ANY", "ANY", "SingleCellExperiment"), function(x, i, j, ..., value) {
    if (missing(i) && missing(j)) {
        int_elementMetadata(x) <- int_elementMetadata(value)
        int_colData(x) <- int_colData(value)
    }

    if (!missing(i)) {
        ii <- .convert_subset_index(i, rownames(x))
        sout <- .standardize_DataFrames(first=int_elementMetadata(x), last=int_elementMetadata(value))
        sout$first[ii,] <- sout$last
        int_elementMetadata(x) <- sout$first
    }

    if (!missing(j)) {
        jj <- .convert_subset_index(j, colnames(x))
        sout <- .standardize_DataFrames(first=int_colData(x), last=int_colData(value), int.col.data=TRUE)
        sout$first[jj,] <- sout$last
        int_colData(x) <- sout$first
    }

    int_metadata(x) <- int_metadata(value)
    callNextMethod()
})
