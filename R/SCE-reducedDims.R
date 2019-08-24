# Getter/setter functions for reducedDims.

.red_key <- "reducedDims"
.unnamed <- "unnamed"

#' @export
#' @importFrom S4Vectors List
#' @importClassesFrom S4Vectors SimpleList
setMethod("reducedDims", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    x <- updateObject(x)
    value <- as(int_colData(x)[[.red_key]], "SimpleList")
    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- colnames(x)
        }
    }
    value
})

#' @export
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame
setReplaceMethod("reducedDims", "SingleCellExperiment", function(x, value) {
    x <- updateObject(x)

    if (length(value)==0L) {
        collected <- int_colData(x)[,0]
    } else {
        nrows <- vapply(value, nrow, FUN.VALUE = 0L)
        if (!all(nrows == ncol(x))) {
            stop("elements of replacement 'reducedDims' do not have the correct number of rows")
        }
        collected <- do.call(DataFrame, lapply(value, I))
        if (is.null(names(value))) {
            colnames(collected) <- paste0(.unnamed, seq_along(value))
        }
    }

    int_colData(x)[[.red_key]] <- collected
    x
})

#' @export
setMethod("reducedDimNames", "SingleCellExperiment", function(x) {
    x <- updateObject(x)
    colnames(int_colData(x)[[.red_key]])
})

#' @export
setReplaceMethod("reducedDimNames", c("SingleCellExperiment", "character"), function(x, value) {
    x <- updateObject(x)
    colnames(int_colData(x)[[.red_key]]) <- value
    x
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE) {

    if (identical(length(reducedDimNames(x)), 0L)) {
        .Deprecated(msg="NULL is deprecated.")
        return(NULL)
        # To deprecate NULL and throw an error instead, remove the two lines above.
        stop(
            "'reducedDim(<", class(x), ">, ...) ",
            "length(reducedDims(<", class(x), ">)) is 0'")
    }

    reducedDim(x, 1L, withDimnames)
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type=1, withDimnames=TRUE) {
    x <- updateObject(x)
    internals <- int_colData(x)[[.red_key]]

    out <- tryCatch({
        internals[, type]
    }, error=function(err) {
        .Deprecated(msg="NULL is deprecated.")
        return(NULL)
        # To deprecate NULL and throw an error instead, remove the two lines above.
        stop("'reducedDim(<", class(x), ">, type=\"numeric\", ...)' ",
             "invalid subscript 'type'\n", conditionMessage(err))
    })

    if (withDimnames) {
        rownames(out) <- colnames(x)
    }

    out
})

#' @export
setMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    msg <- paste0("'reducedDim(<", class(x), ">, type=\"character\", ...)' ",
                  "invalid subscript 'type'")

    x <- updateObject(x)
    internals <- int_colData(x)[[.red_key]]

    out <- tryCatch({
        internals[, type]
    }, error=function(err) {
        .Deprecated(msg="NULL is deprecated.")
        return(NULL)
        # To deprecate NULL and throw an error instead, remove the two lines above.
        stop(msg, "\n'", type, "' not in names(reducedDims(<", class(x), ">))")
    })

    if (withDimnames) {
        rownames(out) <- colnames(x)
    }

    out
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    type <- 1L

    if (0L == length(reducedDimNames(x))){
        ## Implementation 1
        # https://github.com/drisso/SingleCellExperiment/pull/35#issuecomment-522258649
        type <- paste0(.unnamed, 1L)
        ## Implementation 2
        # `reducedDims<-`` above creates character(n) names
        # type <- character(1)
        ## Implementation 3
        # This would implement the same behaviour as SummarizedExperiment::assay<-
        # stop("'reducedDim(<", class(x), ">) <- value' ", "length(reducedDims(<",
        #      class(x), ">)) is 0")
    }

    reducedDim(x, type) <- value
    x
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "numeric"), function(x, type = 1, ..., value) {
    x <- updateObject(x)

    if (!is.null(value) && !identical(nrow(value), ncol(x))) {
        stop("replacement 'reducedDim' has a different number of rows than 'ncol(x)'")
    }

    internals <- int_colData(x)

    if (type[1] > ncol(internals[[.red_key]])) {
        stop("invalid subscript 'type'\nsubscript out of bounds")
    }

    internals[[.red_key]][[type]] <- value
    int_colData(x) <- internals
    x
})

#' @export
setReplaceMethod("reducedDim", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    x <- updateObject(x)

    internals <- int_colData(x)
    if (!is.null(value) && !identical(nrow(value), ncol(x))) {
        stop("replacement 'reducedDim' has a different number of rows than 'ncol(x)'")
    }

    internals[[.red_key]][[type]] <- value
    int_colData(x) <- internals
    x
})
