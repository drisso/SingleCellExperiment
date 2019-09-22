#' @title
#' Named assay getters and setters
#'
#' @description
#' These are methods for getting or setting \code{assay(sce, i=X, ...)} 
#' where \code{sce} is a \linkS4class{SingleCellExperiment} object and \code{X} is the name of the method.
#' For example, \code{counts} will get or set \code{X="counts"}.
#' This provides some convenience for users as well as encouraging standardization of assay names across packages.
#'
#' @section Available methods:
#' In the following code snippets, \code{x} is a \linkS4class{SingleCellExperiment} object,
#' \code{value} is a matrix-like object with the same dimensions as \code{x},
#' and \code{...} are further arguments passed to \code{\link{assay}} (for the getter) or \code{\link{assay<-}} (for the setter).
#' \describe{
#' \item{\code{counts(x, ...)}, \code{counts(x, ...) <- value}:}{
#' Get or set a matrix of raw count data, e.g., number of reads or transcripts.
#' }
#' \item{\code{normcounts(x, ...)}, \code{normcounts(x, ...) <- value}:}{
#' Get or set a matrix of normalized values on the same scale as the original counts.
#' For example, counts divided by cell-specific size factors that are centred at unity.
#' }
#' \item{\code{logcounts(x, ...)}, \code{logcounts(x, ...) <- value}:}{
#' Get or set a matrix of log-transformed counts or count-like values.
#' In most cases, this will be defined as log-transformed \code{normcounts}, e.g., using log base 2 and a pseudo-count of 1.
#' }
#' \item{\code{cpm(x, ...)}, \code{cpm(x, ...) <- value}:}{
#' Get or set a matrix of counts-per-million values.
#' This is the read count for each gene in each cell, divided by the library size of each cell in millions.
#' }
#' \item{\code{tpm(x, ...)}, \code{tpm(x, ...) <- value}:}{
#' Get or set a matrix of transcripts-per-million values.
#' This is the number of transcripts for each gene in each cell, divided by the total number of transcripts in that cell (in millions).
#' }
#' \item{\code{weights(x, ...)}, \code{weights(x, ...) <- value}:}{
#' Get or set a matrix of weights, e.g., observational weights to be used in differential expression analysis.
#' }
#' }
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{assay}} and \code{\link{assay<-}}, for the wrapped methods.
#' 
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#' counts(sce) <- matrix(rnorm(nrow(sce)*ncol(sce)), ncol=ncol(sce))
#' dim(counts(sce))
#' 
#' # One possible way of computing normalized "counts"
#' sf <- 2^rnorm(ncol(sce))
#' sf <- sf/mean(sf)
#' normcounts(sce) <- t(t(counts(sce))/sf)
#' dim(normcounts(sce))
#' 
#' # One possible way of computing log-counts
#' logcounts(sce) <- log2(normcounts(sce)+1)
#' dim(normcounts(sce))
#' 
#' @name SCE-assays
#' @rdname assays
#' @docType methods
#' @aliases
#' counts
#' counts<-
#' counts,SingleCellExperiment-method
#' counts<-,SingleCellExperiment-method
#' normcounts
#' normcounts<-
#' normcounts,SingleCellExperiment-method
#' normcounts<-,SingleCellExperiment-method
#' logcounts
#' logcounts<-
#' logcounts,SingleCellExperiment-method
#' logcounts<-,SingleCellExperiment-method
#' cpm
#' cpm<-
#' cpm,SingleCellExperiment-method
#' cpm<-,SingleCellExperiment-method
#' tpm
#' tpm<-
#' tpm,SingleCellExperiment-method
#' tpm<-,SingleCellExperiment-method
#' weights
#' weights<-
#' weights,SingleCellExperiment-method
#' weights<-,SingleCellExperiment-method
NULL

GET_FUN <- function(exprs_values, ...) {
    (exprs_values) # To ensure evaluation
    function(object, ...) {
        assay(object, i=exprs_values, ...)
    }
}

SET_FUN <- function(exprs_values, ...) {
    (exprs_values) # To ensure evaluation
    function(object, ..., value) {
        assay(object, i=exprs_values, ...) <- value
        object
    }
}

#' @export
#' @importFrom BiocGenerics counts
setMethod("counts", "SingleCellExperiment", GET_FUN("counts"))

#' @export
#' @importFrom BiocGenerics "counts<-"
setReplaceMethod("counts", c("SingleCellExperiment", "ANY"), SET_FUN("counts"))

#' @export
setMethod("logcounts", "SingleCellExperiment", GET_FUN("logcounts"))

#' @export
setReplaceMethod("logcounts", c("SingleCellExperiment", "ANY"), SET_FUN("logcounts"))

#' @export
setMethod("normcounts", "SingleCellExperiment", GET_FUN("normcounts"))

#' @export
setReplaceMethod("normcounts", c("SingleCellExperiment", "ANY"), SET_FUN("normcounts"))

#' @export
setMethod("cpm", "SingleCellExperiment", GET_FUN("cpm"))

#' @export
setReplaceMethod("cpm", c("SingleCellExperiment", "ANY"), SET_FUN("cpm"))

#' @export
setMethod("tpm", "SingleCellExperiment", GET_FUN("tpm"))

#' @export
setReplaceMethod("tpm", c("SingleCellExperiment", "ANY"), SET_FUN("tpm"))

#' @export
#' @importFrom BiocGenerics weights
setMethod("weights", "SingleCellExperiment", GET_FUN("weights"))

#' @export
setReplaceMethod("weights", c("SingleCellExperiment", "ANY"), SET_FUN("weights"))
