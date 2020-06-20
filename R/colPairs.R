#' Column pair methods
#'
#' Methods to get or set column pairings in a \linkS4class{SingleCellExperiment} object.
#' These are typically used to store and retrieve relationships between cells, 
#' e.g., in nearest-neighbor graphs or for inferred cell-cell interactions.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{colPair(x, type, asSparse=FALSE)}:}{
#' Retrieves a \linkS4class{SelfHits} object where each entry represents a pair of columns of \code{x}
#' and has number of nodes equal to \code{ncol(x)}.
#' Associated data values are stored in the \code{"value"} field in \code{\link{mcols}}.
#' \code{type} is either a string specifying the name of the column pairing in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired result.
#' If \code{asSparse=FALSE}, a sparse matrix is returned instead where each row/column corresponds to a column of \code{x}.
#' }
#' \item{\code{colPairNames(x)}:}{
#' Returns a character vector containing the names of all dimensionality reduction results in \code{x}.
#' This is guaranteed to be of the same length as the number of results, though the names may not be unique.
#' }
#' \item{\code{colPairs(x, asSparse=FALSE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more column pairings as \linkS4class{SelfHits} objects.
#' If \code{asSparse=FALSE}, each entry is instead a sparse matrix.
#' }
#' }
#'
#' @section Single setter:
#' \code{colPair(x, type) <- value} will add or replace a column pairing
#' in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{type} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{type} is missing, \code{value} is assigned to the first result.
#' If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed1"}.
#' \item If \code{type} is a numeric scalar, it must be within the range of existing results, and \code{value} will be assigned to the result at that index.
#' \item If \code{type} is a string and a result exists with this name, \code{value} is assigned to to that result.
#' Otherwise a new result with this name is append to the existing list of results.
#' }
#'
#' \code{value} is expected to be a \linkS4class{SelfHits} with number of nodes equal to \code{ncol(x)}.
#' It may also be a sparse matrix with number of rows and columns equal to \code{ncol(x)}.
#' Alternatively, if \code{value} is \code{NULL}, the result corresponding to \code{type} is removed from \code{x}.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{colPairs(x) <- value}:}{
#' Replaces all column pairings in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of \linkS4class{SelfHits} or sparse matrices,
#' each of which is subject to the constraints described for the single setter.
#'
#' If \code{value} is named, those names will be used to name the dimensionality reduction results in \code{x}.
#' Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed"}.
#'
#' If \code{value} is \code{NULL}, all dimensionality reduction results in \code{x} are removed.
#' }
#' \item{\code{colPairNames(x) <- value}:}{
#' Replaces all names for column pairings in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of results currently in \code{x}.
#' }
#' }
#'
#' @section A note on subset replacement:
#' When column-subset replacement is performed on a SingleCellExperiment object (i.e., \code{x[,i] <- y}),
#' pairings in \code{x} are only replaced if both columns belong in \code{i}.
#' Pairings with only one column in \code{i} are preserved in order to ensure that \code{x[,i] <- x[,i]} is a no-op.
#' However, if we are replacing the identity of the features in \code{x[i,]},
#' it is unlikely that the pairings involving the old identities are applicable to the replacement features in \code{y}.
#' In such cases, additional pruning may be required to remove all pairs involving \code{i} prior to replacement.
#'
#' @author Aaron Lun 
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE)
#'
#' # Making up some regulatory pairings:
#' hits <- SelfHits(
#'     sample(ncol(sce), 10),
#'     sample(ncol(sce), 10),
#'     nnode=ncol(sce)
#' )
#' mcols(hits)$value <- runif(10)
#'
#' colPair(sce, "regulators") <- hits
#' colPair(sce, "regulators")
#' 
#' as.mat <- colPair(sce, "regulators", asSparse=TRUE)
#' class(as.mat)
#'
#' colPair(sce, "coexpression") <- hits
#' colPairs(sce)
#'
#' colPair(sce, "regulators") <- NULL
#' colPairs(sce)
#'
#' colPairs(sce) <- SimpleList()
#' colPairs(sce)
#'
#' @seealso
#' \code{\link{rowPairs}}, for the row equivalent.
#'
#' @name colPairs
#' @docType methods
#' @aliases
#' colPair colPairs colPairNames
#' colPair,SingleCellExperiment,missing-method
#' colPair,SingleCellExperiment,numeric-method
#' colPair,SingleCellExperiment,character-method
#' colPairs,SingleCellExperiment-method
#' colPairNames,SingleCellExperiment-method
#' colPair<- colPairs<- colPairNames<-
#' colPair<-,SingleCellExperiment,missing-method
#' colPair<-,SingleCellExperiment,numeric-method
#' colPair<-,SingleCellExperiment,character-method
#' colPairs<-,SingleCellExperiment-method
#' colPairNames<-,SingleCellExperiment,character-method
NULL

.colp_key <- "colPairs"

#' @export
#' @importFrom S4Vectors mcols endoapply
setMethod("colPairs", "SingleCellExperiment", function(x, asSparse=FALSE) {
    value <- .get_internal_all(x, 
        getfun=int_colData, 
        key=.colp_key)

    value <- endoapply(value, .get_hits)
    if (asSparse) {
        for (i in seq_along(value)) {
            value[[i]] <- .hits2mat(value[[i]], mcols(value[[i]])$value)
        }
    }

    value
})

#' @importClassesFrom S4Vectors SelfHits
.any2dualsubset <- function(value) {
    if (!is(value, "SelfHits")) {
        value <- .mat2hits(value)
    }
    DualSubset(value)
}

#' @export
setReplaceMethod("colPairs", "SingleCellExperiment", function(x, value) {
    .set_internal_all(x, value, 
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.colp_key,
        convertfun=.any2dualsubset,
        xdimfun=ncol,
        vdimfun=length,
        funstr="colPairs",
        xdimstr="ncol",
        vdimstr="nodes")
})

#' @export
setMethod("colPairNames", "SingleCellExperiment", function(x) {
    .get_internal_names(x, 
        getfun=int_colData, 
        key=.colp_key)
})

#' @export
setReplaceMethod("colPairNames", c("SingleCellExperiment", "character"), function(x, value) {
    .set_internal_names(x, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.colp_key)
})

#' @export
setMethod("colPair", c("SingleCellExperiment", "missing"), function(x, type, asSparse=FALSE) {
    .get_internal_missing(x, 
        basefun=colPair, 
        namefun=colPairNames, 
        funstr="colPair",
        asSparse=asSparse)
})

#' @export
setMethod("colPair", c("SingleCellExperiment", "numeric"), function(x, type, asSparse=FALSE) {
    out <- .get_internal_integer(x, type,
        getfun=int_colData,
        key=.colp_key,
        funstr="colPair",
        substr="type")

    out <- .get_hits(out)
    if (asSparse) {
        out <- .hits2mat(out)
    }

    out
})

#' @export
setMethod("colPair", c("SingleCellExperiment", "character"), function(x, type, asSparse=FALSE) {
    out <- .get_internal_character(x, type,
        getfun=int_colData,
        key=.colp_key,
        funstr="colPair",
        substr="type",
        namestr="colPairNames")

    out <- .get_hits(out)
    if (asSparse) {
        out <- .hits2mat(out)
    }

    out
})

#' @export
setReplaceMethod("colPair", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    .set_internal_missing(x, value,
        basefun=`colPair<-`,
        namefun=colPairNames
    )
})

#' @export
setReplaceMethod("colPair", c("SingleCellExperiment", "numeric"), function(x, type, ..., value) {
    .set_internal_numeric(x, type, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.colp_key,
        convertfun=.any2dualsubset,
        xdimfun=ncol,
        vdimfun=length,
        funstr="colPair",
        xdimstr="ncol",
        vdimstr="nodes",
        substr="type")
})

#' @export
setReplaceMethod("colPair", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    .set_internal_character(x, type, value, 
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.colp_key,
        convertfun=.any2dualsubset,
        xdimfun=ncol, 
        vdimfun=length,
        funstr="colPair", 
        xdimstr="ncol",
        vdimstr="nodes", 
        substr="type") 
})
