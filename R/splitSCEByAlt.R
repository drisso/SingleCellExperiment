#' @export
#' @importFrom SummarizedExperiment colData colData<-
splitSCEByAlt <- function(x, f, ref=NULL) {
    by.feat <- split(seq_along(f), f)
    if (is.null(ref)) {
        ref <- names(by.feat)[which.max(lengths(by.feat))]
    }

    x0 <- x[by.feat[[ref]],]
    for (other in setdiff(names(by.feat), ref)) {
        # Clearing out the colData() before adding it.
        subset <- x[by.feat[[other]],]
        colData(subset) <- colData(se)[,0]
        altExp(x0, other) <- subset
    }
    x0
}
