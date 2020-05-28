#' Get or set the row subset
#'
#' Get or set the row subset in an instance of a \linkS4class{SingleCellExperiment} class.
#' This is assumed to specify some interesting subset of genes to be favored in downstream analyses.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param field String containing the name of the field in the \code{\link{rowData}} to get or set subsetting data.
#' @param value Any character, logical or numeric vector specifying rows of \code{x} to include in the subset.
#' Alternatively \code{NULL}, in which case existing subsetting information is removed.
#' @param ... Additional arguments, currently ignored.
#' @param onAbsence String indicating an additional action to take when labels are absent:
#' nothing (\code{"none"}), a warning (\code{"warn"}) or an error (\code{"error"}).
#'
#' @details
#' A frequent task in single-cell data analyses is to focus on a subset of genes of interest,
#' e.g., highly variable genes, derived marker genes for clusters, known markers for cell types.
#' A related task is to filter out uninteresting genes such as ribosomal protein genes or mitochondrial transcripts,
#' in which case we want to subset to exclude those genes.
#' 
#' These functions store a set of genes of interest inside a \linkS4class{SingleCellExperiment}
#' for later retrieval and use in downstream functions.
#' Character and numeric \code{value} are converted to logical vectors that are parallel to the rows of \code{x},
#' allowing them to be added to the \code{\link{rowData}} for synchronized row-level operations.
#'
#' For developers, \code{onAbsence} is provided to make it easier to mandate that \code{x} actually has labels.
#' This avoids silent \code{NULL} values that flow to the rest of the function and make debugging difficult.
#'
#' @author Aaron Lun
#'
#' @return
#' For \code{rowSubset}, a logical vector is returned specifying the rows to retain in the subset of interest.
#' If no subset is available, a \code{NULL} is returned (and/or a warning or error, depending on \code{onAbsence}).
#' 
#' For \code{rowSubset<-}, a modified \code{x} is returned a subsetting vector in its \code{\link{rowData}}.
#' 
#' @seealso
#' \linkS4class{SingleCellExperiment}, for the underlying class definition. 
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#'
#' rowSubset(sce, "hvgs") <- 1:10
#' rowSubset(sce, "hvgs")
#' 
#' rowSubset(sce) <- rbinom(nrow(sce), 1, 0.5)==1
#' rowSubset(sce)
#'
#' @docType methods
#' @name rowSubset
#' @aliases rowSubset rowSubset<-
NULL

.subset_field <- "subset"

#' @export
#' @rdname rowSubset
#' @importFrom SummarizedExperiment rowData
setMethod("rowSubset", "SingleCellExperiment", function(x, field="subset", onAbsence="none") {
    output <- rowData(x)[[field]]
    .absent_action(x, val=output, fun="rowSubset", onAbsence=onAbsence)
    output
})

#' @export
#' @rdname rowSubset
#' @importFrom SummarizedExperiment rowData<- rowData
setReplaceMethod("rowSubset", "SingleCellExperiment", function(x, field="subset", ..., value) {
    if (is.numeric(value)) {
        value <- seq_len(nrow(x)) %in% value
    } else if (is.character(value)) {
        value <- rownames(x) %in% value
    }

    rowData(x)[[field]] <- value
    x
})
