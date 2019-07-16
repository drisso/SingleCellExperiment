\name{altExp}

\alias{altExp}
\alias{altExp,SingleCellExperiment-method}
\alias{altExps}
\alias{altExps,SingleCellExperiment-method}
\alias{altExpNames}
\alias{altExpNames,SingleCellExperiment-method}
\alias{altAssays}
\alias{altAssays,SingleCellExperiment-method}
\alias{altAssayNames}
\alias{altAssayNames,SingleCellExperiment-method}
\alias{altAssay}
\alias{altAssay,SingleCellExperiment-method}
\alias{altRowData}
\alias{altRowData,SingleCellExperiment-method}
\alias{altRowNames}
\alias{altRowNames,SingleCellExperiment-method}

\alias{altExp<-}
\alias{altExp<-,SingleCellExperiment-method}
\alias{altExps<-}
\alias{altExps<-,SingleCellExperiment-method}
\alias{altExpNames<-}
\alias{altExpNames<-,SingleCellExperiment-method}
\alias{altAssays<-}
\alias{altAssays<-,SingleCellExperiment-method}
\alias{altAssayNames<-}
\alias{altAssayNames<-,SingleCellExperiment-method}
\alias{altAssay<-}
\alias{altAssay<-,SingleCellExperiment-method}
\alias{altRowData<-}
\alias{altRowData<-,SingleCellExperiment-method}
\alias{altRowNames<-}
\alias{altRowNames<-,SingleCellExperiment-method}

% Dumped here to avoid problems with CHECK.
\alias{[,SummarizedExperimentByColumn,ANY,ANY,ANY-method}
\alias{[<-,SummarizedExperimentByColumn,ANY,ANY,ANY-method}
\alias{c,SummarizedExperimentByColumn-method}
\alias{length,SummarizedExperimentByColumn-method}
\alias{names,SummarizedExperimentByColumn-method}
\alias{names<-,SummarizedExperimentByColumn-method}

\title{Alternative experiments}

\description{Store alternative \linkS4class{SummarizedExperiment} objects inside a \linkS4class{SingleCellExperiment}, containing data for a different set of features but the same samples.
}

\usage{
\S4method{altExps}{SingleCellExperiment}(x)

\S4method{altExpNames}{SingleCellExperiment}(x)

\S4method{altExp}{SingleCellExperiment}(x, e=1)

\S4method{altAssays}{SingleCellExperiment}(x, e=1, ...)

\S4method{altAssayNames}{SingleCellExperiment}(x, e=1, ...)

\S4method{altAssay}{SingleCellExperiment}(x, e=1, i=1, ...)

\S4method{altRowData}{SingleCellExperiment}(x, e=1, ...)

\S4method{altRowNames}{SingleCellExperiment}(x, e=1)

\S4method{altExps}{SingleCellExperiment}(x) <- value

\S4method{altExpNames}{SingleCellExperiment}(x) <- value

\S4method{altExp}{SingleCellExperiment}(x, e=1, ...) <- value

\S4method{altAssays}{SingleCellExperiment}(x, e=1, ...) <- value

\S4method{altAssayNames}{SingleCellExperiment}(x, e=1, ...) <- value

\S4method{altAssay}{SingleCellExperiment}(x, e=1, i=1, ...) <- value

\S4method{altRowData}{SingleCellExperiment}(x, e=1, ...) <- value

\S4method{altRowNames}{SingleCellExperiment}(x, e=1, ...) <- value
}

\arguments{
\item{x}{A \linkS4class{SingleCellExperiment} object.}
\item{e}{An integer scalar or string specifying the alternative experiment to get or set.}
\item{i}{An integer scalar or string specifying the assay of the alternative experiment to get or set.}
\item{...}{Further arguments to be passed to relevant functions.

For \code{altAssays}, these are passed to \code{\link{assays}}.

For \code{altAssayNames}, these are passed to \code{\link{assayNames}}.

For \code{altAssay}, these are passed to \code{\link{assay}}.

For \code{altRowData}, these are passed to \code{\link{rowData}}.

For \code{altAssays<-}, these are passed to \code{\link{assays<-}}.

For \code{altAssayNames<-}, these are passed to \code{\link{assayNames<-}}.

For \code{altAssay<-}, these are passed to \code{\link{assay<-}}.

For \code{altRowData<-}, these are passed to \code{\link{rowData<-}}.

For other methods, these are ignored.
}
\item{value}{An appropriate replacement value.

For \code{altExpNames<-}, a character vector of length equal to the number of alternative experiments.

For \code{altExps<-}, a named list-like object where each element is a \linkS4class{SummarizedExperiment} object with the same number of rows.
Alternatively, a \code{NULL} value to remove all alternativeExperiments from \code{x}.

For \code{altExp<-}, a SummarizedExperiment instance with number of samples equal to \code{ncol(x)}.
Alternatively, a \code{NULL} value to remove the Experiment specified by \code{e}.

For \code{altAssays<-}, a \linkS4class{Assays} object where each element has number of columns equal to \code{ncol(x)}.

For \code{altAssayNames<-}, a character vector of length equal to the number of assays.

For \code{altAssay<-}, a matrix-like object with dimensions equal to the alternative experiment specified by \code{e}.

For \code{altRowData<-}, a \linkS4class{DataFrame} with number of rows equal to the alternative experiment specified by \code{e}.

For \code{altRowNames<-}, a character vector of length equal to the number of rows of the alternative experiment specified by \code{e}.
Alternatively, \code{NULL}.
}
}

\details{
In some experiments, different features must be normalized differently or have different row-level metadata.
Typical examples would be for spike-in transcripts in plate-based experiments and antibody or CRISPR tags in CITE-seq experiments.
These data cannot be stored in the main \code{assays} of the \linkS4class{SingleCellExperiment} itself.
However, it is still desirable to store these features \emph{somewhere} in the SingleCellExperiment.
This simplifies account-keeping in long workflows and ensure that samples remain synchronised.

To facilitate this, the \linkS4class{SingleCellExperiment} allows for \dQuote{alternative Experiments}.
This stores nested \linkS4class{SummarizedExperiment}-class objects inside a SingleCellExperiment object \code{x}, in a manner that guarantees that the nested objects have the same columns in the same order as those in \code{x}.
Methods are provided to enable convenient access to and manipulation of these alternative Experiments.

Note that an \code{altColData} method is \emph{not} provided.
Column-level metadata should rightfully be stored in the \code{colData} of the \linkS4class{SingleCellExperiment}, given that the same set of samples is shared across the main and nested Experiment objects.
}

\value{
\code{altExpNames} returns a character vector of names for the alternative Experiments.

\code{altExps} returns a named \linkS4class{SimpleList} of SummarizedExperiment objects.

\code{altExp} returns a SummarizedExperiment instance or a subclass thereof.

\code{altAssays} returns a \linkS4class{Assays} instance. 

\code{altAssayNames} returns a character vector of assay names.

\code{altAssay} returns a matrix-like object, usually containing experimental observations for each alternative feature.

\code{altRowData} returns a \linkS4class{DataFrame} of row-level metadata. 

\code{altRowNames} returns a character vector of row names or \code{NULL}.

The replacement methods return a SingleCellExperiment with the same dimensions as \code{x} and modified alternative Experiments.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{splitSCEByAlt}}, for a convenient way of adding alternative Experiments from existing features.
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
dim(counts(sce))

# Mocking up some alternative Experiments.
se1 <- SummarizedExperiment(matrix(rpois(1000, 5), ncol=ncol(se)))
rowData(se1)$stuff <- sample(LETTERS, nrow(se1), replace=TRUE)
se2 <- SummarizedExperiment(matrix(rpois(500, 5), ncol=ncol(se)))
rowData(se2)$blah <- sample(letters, nrow(se2), replace=TRUE)

# Setting the alternative Experiments.
altExp(sce, "spike-in") <- se1
altExp(sce, "CRISPR") <- se2
altRowNames(sce, 2) <- paste0("TAG_", seq_len(nrow(se2)))

# Getting alternative Experimental data.
altExpNames(sce)
altExp(sce, "spike-in")
altExp(sce, 2)

altAssay(sce, "CRISPR", 1)
altRowData(sce, "CRISPR")
altRowNames(sce, "CRISPR")

# Setting alternative Experimental data.
altExpNames(sce) <- c("ERCC", "Ab")
altExp(sce, "ERCC") <- se1[1:2,]

altAssay(sce, "Ab", 2) <- matrix(rnorm(500), ncol=ncol(se2))
altRowData(sce, "Ab")$blah <- toupper(altRowData(sce, "Ab")$blah)
altRowNames(sce, "Ab") <- NULL
}
