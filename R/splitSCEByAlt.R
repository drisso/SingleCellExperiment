#' @export
splitSCEByAlt <- function(x, f, ref=NULL) {
    by.feat <- split(seq_along(f), f)
    if (is.null(ref)) {
        ref <- names(by.feat)[which.max(lengths(by.feat))]
    }

    x0 <- x[by.feat[[ref]],]
    for (other in setdiff(names(by.feat), ref)) {
        altExp(x0, other) <- x[by.feat[[other]],]
    }
    x0
}
