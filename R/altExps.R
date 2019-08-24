#' Alternative Experiment methods
#'
#' @description
#' In some experiments, different features must be normalized differently or have different row-level metadata.
#' Typical examples would be for spike-in transcripts in plate-based experiments and antibody or CRISPR tags in CITE-seq experiments.
#' These data cannot be stored in the main \code{assays} of the \linkS4class{SingleCellExperiment} itself.
#' However, it is still desirable to store these features \emph{somewhere} in the SingleCellExperiment.
#' This simplifies book-keeping in long workflows and ensure that samples remain synchronised.
#' 
#' To facilitate this, the \linkS4class{SingleCellExperiment} class allows for \dQuote{alternative Experiments}.
#' Nested \linkS4class{SummarizedExperiment}-class objects are stored inside the SingleCellExperiment object \code{x}, in a manner that guarantees that the nested objects have the same columns in the same order as those in \code{x}.
#' Methods are provided to enable convenient access to and manipulation of these alternative Experiments.
#' Each alternative Experiment should contain experimental data and row metadata for a distinct set of features.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{altExp(x, e, withColData=TRUE)}:}{
#' Retrieves a \linkS4class{SummarizedExperiment} containing alternative features (rows) for all cells (columns) in \code{x}.
#' \code{e} is either a string specifying the name of the alternative Experiment in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired Experiment.
#' If \code{withColData=TRUE}, the column metadata of the output object is set to \code{\link{colData}(x)}.
#' }
#' \item{\code{altExpNames(x)}:}{
#' Returns a character vector containing the names of all alternative Experiments in \code{x}.
#' This is guaranteed to be of the same length as the number of results, though the names may not be unique.
#' }
#' \item{\code{altExps(x, withColData=TRUE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more \linkS4class{SummarizedExperiment} objects.
#' Each object has the same number of columns.
#' If \code{withColData=TRUE}, the column metadata of each output object is set to \code{\link{colData}(x)}.
#' }
#' }
#'
#' @section Single-object setter:
#' \code{altExp(x, e) <- value} will add or replace an alternative Experiment 
#' in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{e} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{e} is missing, \code{value} is assigned to the first result.
#' If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed1"}.
#' \item If \code{e} is a numeric scalar, it must be within the range of existing results, and \code{value} will be assigned to the result at that index.
#' \item If \code{e} is a string and a result exists with this name, \code{value} is assigned to to that result.
#' Otherwise a new result with this name is append to the existing list of results.
#' }
#'
#' \code{value} is expected to be a SummarizedExperiment object with number of columns equal to \code{ncol(x)}.
#' Alternatively, if \code{value} is \code{NULL}, the alternative Experiment at \code{e} is removed from the object.
#' 
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{altExps(x) <- value}:}{
#' Replaces all alterrnative Experiments in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of SummarizedExperiment objects 
#' with number of columns equal to \code{ncol(x)}.
#'
#' If \code{value} is named, those names will be used to name the alternative Experiments in \code{x}.
#' Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed"}.
#' 
#' If \code{value} is \code{NULL}, all alternative Experiments in \code{x} are removed.
#' }
#' \item{\code{altExpNames(x) <- value}:}{
#' Replaces all names for alternative Experiments in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of results currently in \code{x}.
#' }
#' }
#'
#' @seealso
#' \code{\link{splitAltExps}}, for a convenient way of adding alternative Experiments from existing features.
#'
#' @author Aaron Lun
#' 
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#' dim(counts(sce))
#' 
#' # Mocking up some alternative Experiments.
#' se1 <- SummarizedExperiment(matrix(rpois(1000, 5), ncol=ncol(se)))
#' rowData(se1)$stuff <- sample(LETTERS, nrow(se1), replace=TRUE)
#' se2 <- SummarizedExperiment(matrix(rpois(500, 5), ncol=ncol(se)))
#' rowData(se2)$blah <- sample(letters, nrow(se2), replace=TRUE)
#' 
#' # Setting the alternative Experiments.
#' altExp(sce, "spike-in") <- se1
#' altExp(sce, "CRISPR") <- se2
#' 
#' # Getting alternative Experimental data.
#' altExpNames(sce)
#' altExp(sce, "spike-in")
#' altExp(sce, 2)
#' 
#' # Setting alternative Experimental data.
#' altExpNames(sce) <- c("ERCC", "Ab")
#' altExp(sce, "ERCC") <- se1[1:2,]
#' 
#' @name altExps
#' @aliases 
#' altExp altExps altExpNames
#' altExp,SingleCellExperiment,missing-method
#' altExp,SingleCellExperiment,numeric-method
#' altExp,SingleCellExperiment,character-method
#' altExps,SingleCellExperiment-method
#' altExpNames,SingleCellExperiment-method
#' altExp<- altExps<- altExpNames<-
#' altExp<-,SingleCellExperiment,missing-method
#' altExp<-,SingleCellExperiment,numeric-method
#' altExp<-,SingleCellExperiment,character-method
#' altExps<-,SingleCellExperiment-method
#' altExpNames<-,SingleCellExperiment-method
#' [,SummarizedExperimentByColumn,ANY,ANY,ANY-method
#' [<-,SummarizedExperimentByColumn,ANY,ANY,ANY-method
#' c,SummarizedExperimentByColumn-method
#' length,SummarizedExperimentByColumn-method
#' names,SummarizedExperimentByColumn-method
#' names<-,SummarizedExperimentByColumn-method
#'
#' % Dumping the SEBC methods here, so that check doesn't complain.
NULL

.alt_key <- "altExps"

#############################
# Getters.

#' @export
setMethod("altExpNames", "SingleCellExperiment", function(x) {
    colnames(int_colData(x)[[.alt_key]])
})

#' @export
#' @importFrom S4Vectors List
#' @importClassesFrom S4Vectors SimpleList
setMethod("altExps", "SingleCellExperiment", function(x, withColData=FALSE) {
    out <- lapply(int_colData(x)[[.alt_key]], .get_se)
    if (withColData) {
        for (i in seq_along(out)) {
            colData(out[[i]]) <- colData(x)
        }
    }
    as(out, "SimpleList")
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "missing"), function(x, e, withColData=FALSE) {

    if (identical(length(altExpNames(x)), 0L)) {
        stop(
            "'altExp(<", class(x), ">, ...) ",
            "length(altExps(<", class(x), ">)) is 0'")
    }

    altExp(x, 1, withColData)
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "numeric"), function(x, e=1, withColData=FALSE) {
    internals <- int_colData(x)

    out <- tryCatch({
        .get_se(internals[,.alt_key][,e])
    }, error=function(err) {
        stop("'altExp(<", class(x), ">, type=\"numeric\", ...)' ",
             "invalid subscript 'e'\n", conditionMessage(err))
    })

    if (withColData) {
        colData(out) <- colData(x)
    }
    colnames(out) <- colnames(x)
    out
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "character"), function(x, e, withColData=FALSE) {
    msg <- paste0(
        "'altExp(<", class(x), ">, e=\"character\", ...)' ",
        "invalid subscript 'e'")
    internals <- int_colData(x)

    out <- tryCatch({
        .get_se(internals[,.alt_key][,e])
    }, error=function(err) {
        stop(msg, "\n'", e, "' not in names(altExps(<", class(x), ">))")
    })

    if (withColData) {
        colData(out) <- colData(x)
    }
    colnames(out) <- colnames(x)
    out
})

#############################
# Setters.

#' @export
setReplaceMethod("altExpNames", "SingleCellExperiment", function(x, value) {
    colnames(int_colData(x)[[.alt_key]]) <- as.character(value)
    x
})

#' @export
#' @importClassesFrom S4Vectors SimpleList
setReplaceMethod("altExps", "SingleCellExperiment", function(x, value) {
    collected <- int_colData(x)[,0]
    for (i in seq_along(value)) {
        collected[[i]] <- SummarizedExperimentByColumn(value[[i]])
    }
    if (!is.null(names(value))) {
        colnames(collected) <- names(value)
    } else {
        colnames(collected) <- character(length(value))
    }
    int_colData(x)[[.alt_key]] <- collected
    x
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "missing"), function(x, e, ..., value) {
    if (0L == length(altExpNames(x))){
        stop("'altExp(<", class(x), ">) <- value' ", "length(altExps(<",
             class(x), ">)) is 0")
    }
    altExp(x, 1L) <- value
    x
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "numeric"), function(x, e=1, ..., value) {
    internals <- int_colData(x)
    if (e[1] > ncol(internals[[.alt_key]])) {
        stop("invalid subscript 'type'\nsubscript out of bounds")
    }
    if (!is.null(value)) {
        value <- SummarizedExperimentByColumn(value)
    }
    internals[[.alt_key]][[e]] <- value
    int_colData(x) <- internals
    x
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "character"), function(x, e, ..., value) {
    internals <- int_colData(x)
    if (!is.null(value)) {
        value <- SummarizedExperimentByColumn(value)
    }
    internals[[.alt_key]][[e]] <- value
    int_colData(x) <- internals
    x
})
