#' Row pair methods
#'
#' Methods to get or set row pairings in a \linkS4class{SingleCellExperiment} object.
#' These are typically used to store and retrieve relationships between features, 
#' e.g., in gene regulatory or co-expression networks.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{rowPair(x, type, asSparse=FALSE)}:}{
#' Retrieves a \linkS4class{SelfHits} object where each entry represents a pair of rows of \code{x}
#' and has number of nodes equal to \code{nrow(x)}.
#' \code{type} is either a string specifying the name of the row pairing in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired pairing.
#'
#' If \code{asSparse=TRUE}, a sparse matrix is returned instead, see below for details.
#' }
#' \item{\code{rowPairNames(x)}:}{
#' Returns a character vector containing the names of all row pairings in \code{x}.
#' This is guaranteed to be of the same length as the number of pairings, though the names may not be unique.
#' }
#' \item{\code{rowPairs(x, asSparse=FALSE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more row pairings as \linkS4class{SelfHits} objects.
#' If \code{asSparse=FALSE}, each entry is instead a sparse matrix.
#' }
#' }
#'
#' When \code{asSparse=TRUE}, the return value will be a triplet-form sparse matrix 
#' where each row/column corresponds to a row of \code{x}.
#' The values in the matrix will be taken from the first metadata field of the underlying \linkS4class{SelfHits} object,
#' with an error being raised if the first metadata field is not of an acceptable type.
#' If there are duplicate pairs, only the value from the last pair is used. 
#' If no metadata is available, the matrix values are set to \code{TRUE} for all pairs.
#'
#' @section Single setter:
#' \code{rowPair(x, type) <- value} will add or replace a row pairing
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
#' \code{value} is expected to be a \linkS4class{SelfHits} with number of nodes equal to \code{nrow(x)}.
#' Any number of additional fields can be placed in \code{\link{mcols}(value)}.
#' Duplicate row pairs are supported and will not be collapsed into a single entry.
#' 
#' \code{value} may also be a sparse matrix with number of rows and columns equal to \code{nrow(x)}.
#' This is converted into a \linkS4class{SelfHits} object with values stored in the metadata as the \code{"x"} field.
#'
#' Alternatively, if \code{value} is \code{NULL}, the pairings corresponding to \code{type} are removed from \code{x}.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{rowPairs(x) <- value}:}{
#' Replaces all row pairings in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of \linkS4class{SelfHits} or sparse matrices,
#' each of which is subject to the constraints described for the single setter.
#'
#' If \code{value} is named, those names will be used to name the row pairings in \code{x}.
#' Otherwise, unnamed pairings are assigned default names prefixed with \code{"unnamed"}.
#'
#' If \code{value} is \code{NULL}, all row pairings in \code{x} are removed.
#' }
#' \item{\code{rowPairNames(x) <- value}:}{
#' Replaces all names for row pairings in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of pairings currently in \code{x}.
#' }
#' }
#'
#' @section Interaction with SingleCellExperiment operations:
#' When row-subset replacement is performed on a SingleCellExperiment object (i.e., \code{x[i,] <- y}),
#' a pair of rows in \code{rowPair(x)} is only replaced if both rows are present in \code{i}.
#' This replacement not only affects the \code{value} of the pair but also whether it even exists in \code{y}.
#' For example, if a pair exists between two rows in \code{x[i,]} but not in the corresponding rows of \code{y},
#' it is removed upon subset replacement.
#'
#' Importantly, pairs in \code{x} with only one row in \code{i} are preserved by replacement.
#' This ensures that \code{x[i,] <- x[i,]} is a no-op.
#' However, if the replacement is fundamentally altering the identity of the features in \code{x[i,]},
#' it is unlikely that the pairings involving the old identities are applicable to the replacement features in \code{y}.
#' In such cases, additional pruning may be required to remove all pairs involving \code{i} prior to replacement.
#'
#' Another interesting note is that, for some \code{i <- 1:n} where \code{n} is in \code{[1, nrow(x))},
#' \code{rbind(x[i,], x[-i,])} will not return a SingleCellExperiment equal to \code{x} with respect to \code{\link{rowPairs}}.
#' This operation will remove any pairs involving one row in \code{i} and another row outside of \code{i},
#' simply because each individual subset operation will remove pairs involving rows outside of the subset.
#'
#' @author Aaron Lun 
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE)
#'
#' # Making up some regulatory pairings:
#' hits <- SelfHits(
#'     sample(nrow(sce), 10),
#'     sample(nrow(sce), 10),
#'     nnode=nrow(sce)
#' )
#' mcols(hits)$value <- runif(10)
#'
#' rowPair(sce, "regulators") <- hits
#' rowPair(sce, "regulators")
#' 
#' as.mat <- rowPair(sce, "regulators", asSparse=TRUE)
#' class(as.mat)
#'
#' rowPair(sce, "coexpression") <- hits
#' rowPairs(sce)
#'
#' rowPair(sce, "regulators") <- NULL
#' rowPairs(sce)
#'
#' rowPairs(sce) <- SimpleList()
#' rowPairs(sce)
#'
#' @seealso
#' \code{\link{colPairs}}, for the column equivalent.
#'
#' @name rowPairs
#' @docType methods
#' @aliases
#' rowPair rowPairs rowPairNames
#' rowPair,SingleCellExperiment,missing-method
#' rowPair,SingleCellExperiment,numeric-method
#' rowPair,SingleCellExperiment,character-method
#' rowPairs,SingleCellExperiment-method
#' rowPairNames,SingleCellExperiment-method
#' rowPair<- rowPairs<- rowPairNames<-
#' rowPair<-,SingleCellExperiment,missing-method
#' rowPair<-,SingleCellExperiment,numeric-method
#' rowPair<-,SingleCellExperiment,character-method
#' rowPairs<-,SingleCellExperiment-method
#' rowPairNames<-,SingleCellExperiment,character-method
#' [,DualSubset,ANY,ANY,ANY-method
#' [<-,DualSubset,ANY,ANY,ANY-method
#' c,DualSubset-method
#' length,DualSubset-method
NULL

.rowp_key <- "rowPairs"

#' @export
#' @importFrom S4Vectors SimpleList
setMethod("rowPairs", "SingleCellExperiment", function(x, asSparse=FALSE) {
    value <- .get_internal_all(x, 
        getfun=int_elementMetadata, 
        key=.rowp_key)

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
setReplaceMethod("rowPairs", "SingleCellExperiment", function(x, value) {
    .set_internal_all(x, value, 
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.rowp_key,
        convertfun=.any2dualsubset,
        xdimfun=nrow,
        vdimfun=length,
        funstr="rowPairs",
        xdimstr="nrow",
        vdimstr="nodes")
})

#' @export
setMethod("rowPairNames", "SingleCellExperiment", function(x) {
    .get_internal_names(x, 
        getfun=int_elementMetadata, 
        key=.rowp_key)
})

#' @export
setReplaceMethod("rowPairNames", c("SingleCellExperiment", "character"), function(x, value) {
    .set_internal_names(x, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.rowp_key)
})

#' @export
setMethod("rowPair", c("SingleCellExperiment", "missing"), function(x, type, asSparse=FALSE) {
    .get_internal_missing(x, 
        basefun=rowPair, 
        namefun=rowPairNames, 
        funstr="rowPair",
        asSparse=asSparse)
})

#' @export
setMethod("rowPair", c("SingleCellExperiment", "numeric"), function(x, type, asSparse=FALSE) {
    out <- .get_internal_integer(x, type,
        getfun=int_elementMetadata,
        key=.rowp_key,
        funstr="rowPair",
        substr="type")

    out <- .get_hits(out)
    if (asSparse) {
        out <- .hits2mat(out)
    }

    out
})

#' @export
setMethod("rowPair", c("SingleCellExperiment", "character"), function(x, type, asSparse=FALSE) {
    out <- .get_internal_character(x, type,
        getfun=int_elementMetadata,
        key=.rowp_key,
        funstr="rowPair",
        substr="type",
        namestr="rowPairNames")

    out <- .get_hits(out)
    if (asSparse) {
        out <- .hits2mat(out)
    }

    out
})

#' @export
setReplaceMethod("rowPair", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    .set_internal_missing(x, value,
        basefun=`rowPair<-`,
        namefun=rowPairNames
    )
})

#' @export
setReplaceMethod("rowPair", c("SingleCellExperiment", "numeric"), function(x, type, ..., value) {
    .set_internal_numeric(x, type, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.rowp_key,
        convertfun=.any2dualsubset,
        xdimfun=nrow,
        vdimfun=length,
        funstr="rowPair",
        xdimstr="nrow",
        vdimstr="nodes",
        substr="type")
})

#' @export
setReplaceMethod("rowPair", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    .set_internal_character(x, type, value, 
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.rowp_key,
        convertfun=.any2dualsubset,
        xdimfun=nrow, 
        vdimfun=length,
        funstr="rowPair", 
        xdimstr="nrow",
        vdimstr="nodes", 
        substr="type") 
})
