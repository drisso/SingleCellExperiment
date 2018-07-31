.convert_subset_index <- function(subset, names) {
    if (is.character(subset)) {
        fmt <- "index out of bounds: %s"
        subset <- SummarizedExperiment:::.SummarizedExperiment.charbound(subset, names, fmt)
    }
    return(as.vector(subset))
}

#' @importFrom utils head
.standardize_DataFrames <- function(...) {
    all.d <- list(...)
    all.fields <- Reduce(union, lapply(all.d, colnames))

    for (d in seq_along(all.d)) {
        cur.d <- all.d[[d]]
        missing.fields <- setdiff(all.fields, colnames(cur.d))
        if (length(missing.fields)) { # do not forgive; people should fix this.
            lost <- paste0(paste0("'", head(missing.fields, 3), "'"), collapse=", ")
            stop("DataFrame ", d, " does not have ", lost)
        }
        all.d[[d]] <- cur.d[,all.fields,drop=FALSE]
    }

    return(all.d)
}

.standardize_reducedDims <- function(...) {
    args <- list(...)
    all.ncells <- lapply(args, ncol)
    all.rd <- lapply(args, reducedDims)
    all.modes <- Reduce(union, lapply(all.rd, names))

    for (m in all.modes) {
        all.dims <- integer(length(all.rd))
        for (d in seq_along(all.rd)) {
            current <- all.rd[[d]][[m]]
            if (is.null(current)) {
                stop("object ", d, " does not have '", m, "' in 'reducedDims'")
            }
            all.dims[d] <- ncol(current)
        }

        # Checking consistency of dimensions between objects.
        udim <- unique(all.dims)
        if (length(udim)!=1) {
            stop("dimensions of '", m, "' are not consistent between objects")
        }
    }

    # Standardizing the order.
    for (d in seq_along(all.rd)) {
        all.rd[[d]] <- all.rd[[d]][all.modes]
    }
    return(all.rd)
}

