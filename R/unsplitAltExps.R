#' Unsplit the alternative experiments
#'
#' Combine the main and alternative experiments back into one \linkS4class{SingleCellExperiment} object.
#' This is effectively the reverse operation to \code{\link{splitAltExps}}.
#'
#' @param sce A \linkS4class{SingleCellExperiment} containing alternative experiments in the \code{\link{altExps}} slot.
#' @param prefix.rows Logical scalar indicating whether the (non-\code{NULL}) row names should be prefixed with the name of the alternative experiment.
#' @param prefix.cols Logical scalar indicating whether the names of column-related fields should be prefixed with the name of the alternative experiment.
#' If \code{NA}, any \code{\link{colData}} of the \code{\link{altExps}} are ignored.
#' @param delayed Logical scalar indicating whether the combining of the assays should be delayed. 
#'
#' @return A SingleCellExperiment where all features in the alternative experiments of \code{sce} are now features in the main experiment.
#' The output object has no alternative experiments of its own.
#' 
#' @author Aaron Lun
#'
#' @details
#' This function is intended for downstream applications that accept a SingleCellExperiment but are not aware of the \code{\link{altExps}} concept.
#' By consolidating all data together, applications like \pkg{iSEE} can use the same machinery to visualize any feature of interest across all modalities.
#' However, for quantitative analyses, it is usually preferable to keep different modalities separate.
#'
#' Assays with the same name are \code{\link{rbind}}ed together in the output object.
#' If a particular name is not present for any experiment, its values are filled in with the appropriately typed \code{NA} instead.
#' By default, this is done efficiently via \linkS4class{ConstantMatrix} abstractions to avoid actually creating a dense matrix of \code{NA}s.
#' If \code{delayed=FALSE}, the combining of matrices is done without any \linkS4class{DelayedArray} wrappers,
#' yielding a simpler matrix representation at the cost of increasing memory usage.
#'
#' Any \code{\link{colData}} or \code{\link{reducedDims}} in the alternative experiments are added to those of the main experiment.
#' The names of these migrated fields are prefixed by the name of the alternative experiment if \code{prefix.cols=TRUE}.
#' 
#' Setting \code{prefix.rows=FALSE}, \code{prefix.cols=NA} and \code{delayed=FALSE} will reverse the effects of \code{\link{splitAltExps}}.
#'
#' @examples
#' counts <- matrix(rpois(10000, 5), ncol=100)
#' sce <- SingleCellExperiment(assays=list(counts=counts))
#' feat.type <- sample(c("endog", "ERCC", "adt"), nrow(sce),
#'     replace=TRUE, p=c(0.8, 0.1, 0.1))
#' sce <- splitAltExps(sce, feat.type)
#'
#' # Making life a little more complicated.
#' logcounts(sce) <- log2(counts(sce) + 1)
#' sce$cluster <- sample(5, ncol(sce), replace=TRUE)
#' reducedDim(sce, "PCA") <- matrix(rnorm(ncol(sce)*2), ncol=2)
#'  
#' # Now, putting Humpty Dumpty back together again. 
#' restored <- unsplitAltExps(sce)
#' restored
#' 
#' @seealso
#' \code{\link{splitAltExps}}, which does the reverse operation of this function.
#'
#' @export
unsplitAltExps <- function(sce, prefix.rows=TRUE, prefix.cols=TRUE, delayed=TRUE) {
    all.se <- c(list(sce), as.list(altExps(sce)))

    args <- list(
        assays=.unsplit_assays(all.se, delayed=delayed),
        colData=.combine_coldata(all.se, prefix=prefix.cols),
        reducedDims=.combine_reddims(all.se, prefix=prefix.cols),
        metadata=metadata(sce)
    )

    all.rr <- .combine_rowranges(all.se, prefix=prefix.rows)
    if (is(all.rr, "DFrame")) {
        args$rowData <- all.rr
    } else {
        args$rowRanges <- all.rr
    }

    do.call(SingleCellExperiment, args)
}

#' @importFrom DelayedArray DelayedArray ConstantArray
#' @importFrom BiocGenerics rbind
#' @importFrom SummarizedExperiment assayNames assay
.unsplit_assays <- function(all.se, delayed) {
    all_assay_names <- unique(unlist(lapply(all.se, assayNames)))
    combined <- vector("list", length(all_assay_names))
    names(combined) <- all_assay_names

    for (assay_name in all_assay_names) {
        current <- vector("list", length(all.se))

        for (s in seq_along(all.se)) {
            cur.se <- all.se[[s]]
            if (assay_name %in% assayNames(cur.se)) {
                mat <- assay(cur.se, assay_name, withDimnames=FALSE)
                if (delayed) {
                    mat <- DelayedArray(mat)
                }
                current[[s]] <- mat
            } else {
                mat <- ConstantArray(dim(cur.se), value=NA)
                if (!delayed) {
                    mat <- as.matrix(mat)
                }
                current[[s]] <- mat
            }
        }
        a <- do.call(rbind, current)
        rownames(a) <- NULL
        combined[[assay_name]] <- a
    }

    combined
}

#' @importFrom S4Vectors mcols<-
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom SummarizedExperiment rowRanges rowData
.combine_rowranges <- function(all.se, prefix) {
    final.rd <- vector("list", length(all.se))
    for (s in seq_along(all.se)) {
        cur.se <- all.se[[s]]
        if (is(cur.se, "RangedSummarizedExperiment")) {
            final.rd[[s]] <- rowRanges(cur.se)
        }
    }

    # Filling in the empties. We promote everyone to a GRL if any of the individuals are GRLs.
    has.grl <- any(vapply(final.rd, FUN=is, class2="GRangesList", FUN.VALUE=TRUE))
    for (e in which(vapply(final.rd, is.null, FUN.VALUE=TRUE))) {
        cur.se <- all.se[[e]]
        if (has.grl) {
            empty <- GRangesList(rep(list(GRanges()), nrow(cur.se)))
        } else {
            empty <- GRanges(rep("unknown:1-0", nrow(cur.se)))
        }
        mcols(empty) <- rowData(cur.se)
        names(empty) <- rownames(cur.se)
        final.rd[[e]] <- empty
    }

    output <- NULL
    tryCatch(
        output <- do.call(c, final.rd),
        error=function(e) {
            warning(conditionMessage(e))
        }
    )

    if (is.null(output)) {
        output <- lapply(all.se, function(x) rowData(x)[,0])
        output <- do.call(rbind, output)
    }

    if (prefix) {
        prefixes <- rep(c("", sprintf("%s.", names(all.se)[-1])), vapply(all.se, nrow, 0L))
        if (is(output, "DFrame")) {
            if (!is.null(rownames(output))) {
                rownames(output) <- sprintf("%s%s", prefixes, rownames(output))
            }
        } else {
            if (!is.null(names(output))) {
                names(output) <- sprintf("%s%s", prefixes, names(output))
            }
        }
    }

    output
}

#' @importFrom SummarizedExperiment colData
.combine_coldata <- function(all.se, prefix) {
    if (!is.na(prefix)) {
        final.cd <- lapply(all.se, colData)
        final.cd <- unname(final.cd)
        if (prefix) {
            for (i in seq_along(all.se)[-1]) {
                colnames(final.cd[[i]]) <- sprintf("%s.%s", names(all.se)[i], colnames(final.cd[[i]]))
            }
        }
        do.call(cbind, final.cd)
    } else {
        colData(all.se[[1]])
    }
}

.combine_reddims <- function(all.se, prefix) {
    if (!is.na(prefix)) {
        final.rd <- lapply(all.se, function(x) {
            if (is(x, "SingleCellExperiment")) {
                reducedDims(x)
            } else {
                list()
            }
        })

        if (prefix) {
            for (i in seq_along(all.se)[-1]) {
                names(final.rd[[i]]) <- sprintf("%s.%s", names(all.se)[i], names(final.rd[[i]]))
            }
        }
        do.call(c, final.rd)
    } else {
        reducedDims(all.se[[1]])
    }
}
