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
#' \code{type} is either a string specifying the name of the column pairing in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired result.
#'
#' If \code{asSparse=TRUE}, a sparse matrix is returned instead, see below for details.
#' }
#' \item{\code{colPairNames(x)}:}{
#' Returns a character vector containing the names of all column pairings in \code{x}.
#' This is guaranteed to be of the same length as the number of results, though the names may not be unique.
#' }
#' \item{\code{colPairs(x, asSparse=FALSE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more column pairings as \linkS4class{SelfHits} objects.
#' If \code{asSparse=FALSE}, each entry is instead a sparse matrix.
#' }
#' }
#'
#' When \code{asSparse=TRUE}, the return value will be a triplet-form sparse matrix 
#' where each row/column corresponds to a column of \code{x}.
#' The values in the matrix will be taken from the first metadata field of the underlying \linkS4class{SelfHits} object,
#' with an error being raised if the first metadata field is not of an acceptable type.
#' If there are duplicate pairs, only the value from the last pair is used. 
#' If no metadata is available, the matrix values are set to \code{TRUE} for all pairs.
#'
#' @section Single setter:
#' \code{colPair(x, type) <- value} will add or replace a column pairing
#' in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{type} determines how the pairing is added or replaced:
#' \itemize{
#' \item If \code{type} is missing, \code{value} is assigned to the first pairing.
#' If the pairing already exists, its name is preserved; otherwise it is given a default name \code{"unnamed1"}.
#' \item If \code{type} is a numeric scalar, it must be within the range of existing pairings, 
#' and \code{value} will be assigned to the pairing at that index.
#' \item If \code{type} is a string and a pairing exists with this name, \code{value} is assigned to to that pairing.
#' Otherwise a new pairing with this name is append to the existing list of pairings.
#' }
#'
#' \code{value} is expected to be a \linkS4class{SelfHits} with number of nodes equal to \code{ncol(x)}.
#' Any number of additional fields can be placed in \code{\link{mcols}(value)}.
#' Duplicate column pairs are supported and will not be collapsed into a single entry.
#' 
#' \code{value} may also be a sparse matrix with number of rows and columns equal to \code{ncol(x)}.
#' This is converted into a \linkS4class{SelfHits} object with values stored in the metadata as the \code{"x"} field.
#' 
#' Alternatively, if \code{value} is \code{NULL}, the pairings corresponding to \code{type} are removed from \code{x}.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{colPairs(x) <- value}:}{
#' Replaces all column pairings in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of \linkS4class{SelfHits} or sparse matrices,
#' each of which is subject to the constraints described for the single setter.
#'
#' If \code{value} is named, those names will be used to name the column pairings in \code{x}.
#' Otherwise, unnamed pairings are assigned default names prefixed with \code{"unnamed"}.
#'
#' If \code{value} is \code{NULL}, all column pairings in \code{x} are removed.
#' }
#' \item{\code{colPairNames(x) <- value}:}{
#' Replaces all names for column pairings in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of pairings currently in \code{x}.
#' }
#' }
#'
#' @section Interaction with SingleCellExperiment operations:
#' When column-subset replacement is performed on a SingleCellExperiment object (i.e., \code{x[,i] <- y}),
#' a pair of columns in \code{colPair(x)} is only replaced if both columns are present in \code{i}.
#' This replacement not only affects the \code{value} of the pair but also whether it even exists in \code{y}.
#' For example, if a pair exists between two columns in \code{x[,i]} but not in the corresponding columns of \code{y},
#' it is removed upon subset replacement.
#'
#' Importantly, pairs in \code{x} with only one column in \code{i} are preserved by replacement.
#' This ensures that \code{x[,i] <- x[,i]} is a no-op.
#' However, if the replacement is fundamentally altering the identity of the features in \code{x[,i]},
#' it is unlikely that the pairings involving the old identities are applicable to the replacement features in \code{y}.
#' In such cases, additional pruning may be required to remove all pairs involving \code{i} prior to replacement.
#'
#' Another interesting note is that, for some \code{i <- 1:n} where \code{n} is in \code{[1, ncol(x))},
#' \code{cbind(x[,i], x[,-i])} will not return a SingleCellExperiment equal to \code{x} with respect to \code{\link{colPairs}}.
#' This operation will remove any pairs involving one column in \code{i} and another column outside of \code{i},
#' simply because each individual subset operation will remove pairs involving columns outside of the subset.
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
#' @importFrom S4Vectors SimpleList
setMethod("colPairs", "SingleCellExperiment", function(x, asSparse=FALSE) {
    value <- .get_internal_all(x, 
        getfun=int_colData, 
        key=.colp_key)

    value <- lapply(value, .get_hits)
    if (asSparse) {
        for (i in seq_along(value)) {
            value[[i]] <- .hits2mat(value[[i]])
        }
    }

    SimpleList(value)
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
