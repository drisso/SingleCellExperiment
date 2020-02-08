.convert_subset_index <- function(subset, names) {
    if (is.character(subset)) {
        fmt <- "index out of bounds: %s"
        subset <- SummarizedExperiment:::.SummarizedExperiment.charbound(subset, names, fmt)
    }
    return(as.vector(subset))
}

.unnamed <- "unnamed"

.absent_action <- function(object, val, fun, onAbsence=c("none", "warn", "error")) {
    if (is.null(val)) {
        onAbsence <- match.arg(onAbsence)
        if (onAbsence!="none") {
            msg <- sprintf("'%s(<%s>)' returns 'NULL'", fun, class(object)[1])
            if (onAbsence=="warn") {
                warning(msg)
            } else {
                stop(msg)
            }
        }
    }
}
