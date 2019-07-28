.convert_subset_index <- function(subset, names) {
    if (is.character(subset)) {
        fmt <- "index out of bounds: %s"
        subset <- SummarizedExperiment:::.SummarizedExperiment.charbound(subset, names, fmt)
    }
    return(as.vector(subset))
}

#' @importFrom utils head
.standardize_DataFrames <- function(..., int.col.data=FALSE) {
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

    if (int.col.data) {
        if (.red_key %in% all.fields) {
            all.rd <- lapply(all.d, "[[", i=.red_key)
            all.rd <- .standardize_reducedDims(all.rd)
            for (d in seq_along(all.d)) {
                all.d[[d]][[.red_key]] <- all.rd[[d]]
            }
        }
        if (.alt_key %in% all.fields) {
            all.alt <- lapply(all.d, "[[", i=.alt_key)
            all.alt <- .standardize_altExps(all.alt)
            for (d in seq_along(all.d)) {
                all.d[[d]][[.alt_key]] <- all.alt[[d]]
            }
        }
    }

    all.d
}

.standardize_reducedDims <- function(all.rd) {
    all.modes <- Reduce(union, lapply(all.rd, colnames))

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
    all.rd
}

.standardize_altExps <- function(all.alt) {
    all.modes <- Reduce(union, lapply(all.alt, colnames))

    for (m in all.modes) {
        all.dims <- integer(length(all.alt))
        for (d in seq_along(all.alt)) {
            current <- all.alt[[d]][[m]]
            if (is.null(current)) {
                stop("object ", d, " does not have '", m, "' in 'altExps'")
            }
            all.dims[d] <- nrow(.get_se(current))
        }

        # Checking consistency of dimensions between objects.
        udim <- unique(all.dims)
        if (length(udim)!=1) {
            stop("dimensions of '", m, "' are not consistent between objects")
        }
    }

    # Standardizing the order.
    for (d in seq_along(all.alt)) {
        all.alt[[d]] <- all.alt[[d]][all.modes]
    }
    all.alt
}
