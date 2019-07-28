#' @export
#' @importFrom BiocGenerics rbind cbind
setMethod("cbind", "SingleCellExperiment", function(..., deparse.level=1) {
    args <- list(...)
    args <- lapply(args, updateObject)
    int_m <- do.call(c, lapply(args, int_metadata))

    tryCatch({
        int_cd <- do.call(rbind, lapply(args, .filled_int_colData))
    }, error=function(err) {
        stop(paste0("failed to combine 'int_colData'\n", err))
    })

    # Creating a shell to avoid having to pull out .cbind.DataFrame 
    # to fuse metadata along the dimension not being combined.
    row_shells <- lapply(args, .create_shell_rowdata)
    tryCatch({
        combined <- do.call(cbind, row_shells)
    }, error=function(err) {
        stop(paste0("failed to combine 'int_elementMetadata'\n", err))
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
    int_m <- do.call(c, lapply(args, int_metadata))

    tryCatch({
        int_em <- do.call(rbind, lapply(args, int_elementMetadata))
    }, error=function(err) {
        stop(paste0("failed to combine 'int_elementMetadata'\n", err))
    })

    # Creating a shell to avoid having to pull out .cbind.DataFrame 
    # to fuse metadata along the dimension not being combined.
    col_shells <- lapply(args, .create_shell_coldata)
    tryCatch({
        combined <- do.call(rbind, col_shells)
    }, error=function(err) {
        stop(paste0("failed to combine 'int_colData'\n", err))
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
    df <- .filled_int_colData(x)
    SummarizedExperiment(colData=df)
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.create_shell_rowdata <- function(x) {
    SummarizedExperiment(rowData=int_elementMetadata(x))
}

