# Size factor getter/setter functions.

.get_sf_field <- function(type) {
    search <- "size_factor"
    if (!is.null(type)) {
        .Deprecated(msg="'type=' is deprecated.")
        if (length(type)!=1L) {
            stop("'type' must be a character vector of length 1")
        }
        search <- paste0(search, "_", type)
    }
    return(search)
}

#' @export
#' @importFrom BiocGenerics sizeFactors
setMethod("sizeFactors", "SingleCellExperiment", function(object, type=NULL) {
    field <- .get_sf_field(type)
    return(int_colData(object)[[field]])
})

#' @export
#' @importFrom BiocGenerics "sizeFactors<-"
setReplaceMethod(
    f = "sizeFactors",
    signature = signature(
        object = "SingleCellExperiment",
        value = "numeric"
    ),
    signature = function(object, type=NULL, ..., value) {
        field <- .get_sf_field(type)
        cd <- int_colData(object)
        cd[[field]] <- value
        int_colData(object) <- cd

        if (!is.null(type)) {
            .Deprecated(msg="'type=' is deprecated.")
            md <- int_metadata(object)
            if (is.null(value)) {
                md$size_factor_names <- setdiff(md$size_factor_names, type)
            } else {
                md$size_factor_names <- union(md$size_factor_names, type)
            }
            int_metadata(object) <- md
        }
        return(object)
    }
)

#' @export
setMethod("clearSizeFactors", "SingleCellExperiment", function(object) {
    .Deprecated()
    sizeFactors(object) <- NULL

    cd <- int_colData(object)
    for (sf in sizeFactorNames(object)) {
        field <- .get_sf_field(sf)
        cd[[field]] <- NULL
    }
    int_colData(object) <- cd

    md <- int_metadata(object)
    md$size_factor_names <- character(0)
    int_metadata(object) <- md
    return(object)
})

#' @export
setMethod("sizeFactorNames", "SingleCellExperiment", function(object) {
    .Deprecated()
    int_metadata(object)$size_factor_names
})
