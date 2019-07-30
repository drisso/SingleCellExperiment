#' @export
setMethod("objectVersion", "SingleCellExperiment", function(x) {
    int_metadata(x)$version
})

# Sets the show method.

scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

.sce_show <- function(object) {
    callNextMethod()
    scat("reducedDimNames(%d): %s\n", reducedDimNames(object))
    scat("spikeNames(%d): %s\n", suppressWarnings(spikeNames(object)))
    scat("altExpNames(%d): %s\n", altExpNames(object))
}

#' @export
setMethod("show", "SingleCellExperiment", .sce_show)
