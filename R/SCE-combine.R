#' @export
#' @importFrom S4Vectors metadata
#' @importFrom BiocGenerics cbind
setMethod("cbind", "SingleCellExperiment", function(..., deparse.level=1) {
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    out <- callNextMethod()

    internals <- .internal_combiner(list(...), FUN=cbind)
    int_cd <- internals$colData
    int_em <- internals$elementMetadata
    int_m <- internals$metadata

    initialize(out, int_colData=int_cd, int_elementMetadata=int_em, int_metadata=int_m)
})

#' @export
#' @importFrom BiocGenerics rbind
setMethod("rbind", "SingleCellExperiment", function(..., deparse.level=1) {
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }
    out <- callNextMethod()

    internals <- .internal_combiner(list(...), FUN=rbind)
    int_cd <- internals$colData
    int_em <- internals$elementMetadata
    int_m <- internals$metadata

    initialize(out, int_colData=int_cd, int_elementMetadata=int_em, int_metadata=int_m)
})

#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
.internal_combiner <- function(args, FUN) 
# Fixing internal fields. To ensure consistent behaviour with SE's c/rbind,
# we create 'shell' SEs where the internal fields are set to the visible 
# row/colData, c/rbind the SEs, and extract the row/colData back out.
{
    meta_shells <- lapply(args, .create_shell_metadata)
    tryCatch({
        combined <- do.call(FUN, meta_shells)
    }, error=function(err) {
        stop(paste0("failed to combine 'int_metadata'\n", err))
    })
    int_m <- metadata(combined)

    col_shells <- lapply(args, .create_shell_coldata)
    tryCatch({
        combined <- do.call(FUN, col_shells)
    }, error=function(err) {
        stop(paste0("failed to combine 'int_colData'\n", err))
    })
    int_cd <- colData(combined)

    row_shells <- lapply(args, .create_shell_rowdata)
    tryCatch({
        combined <- do.call(FUN, row_shells)
    }, error=function(err) {
        stop(paste0("failed to combine 'int_elementMetadata'\n", err))
    })
    int_em <- rowData(combined)

    list(colData=int_cd, elementMetadata=int_em, metadata=int_m)
}
