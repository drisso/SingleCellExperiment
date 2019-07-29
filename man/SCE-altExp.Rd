\name{altExp}

\alias{altExp}
\alias{altExp,SingleCellExperiment-method}
\alias{altExps}
\alias{altExps,SingleCellExperiment-method}
\alias{altExpNames}
\alias{altExpNames,SingleCellExperiment-method}

\alias{altExp<-}
\alias{altExp<-,SingleCellExperiment-method}
\alias{altExps<-}
\alias{altExps<-,SingleCellExperiment-method}
\alias{altExpNames<-}
\alias{altExpNames<-,SingleCellExperiment-method}

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
\S4method{altExps}{SingleCellExperiment}(x, withColData=TRUE)

\S4method{altExpNames}{SingleCellExperiment}(x)

\S4method{altExp}{SingleCellExperiment}(x, e=1, withColData=TRUE)

\S4method{altExps}{SingleCellExperiment}(x) <- value

\S4method{altExpNames}{SingleCellExperiment}(x) <- value

\S4method{altExp}{SingleCellExperiment}(x, e=1, ...) <- value
}

\arguments{
\item{x}{A \linkS4class{SingleCellExperiment} object.}
\item{withColData}{Logical scalar indicating whether the column metadata of \code{x} should replace that of the output object.}
\item{e}{An integer scalar or string specifying the alternative experiment to get or set.}
\item{...}{Ignored.}
\item{value}{An appropriate replacement value.

For \code{altExpNames<-}, a character vector of length equal to the number of alternative experiments.

For \code{altExps<-}, a named list-like object where each element is a \linkS4class{SummarizedExperiment} object with the same number of columns.
Alternatively, a \code{NULL} value to remove all alternativeExperiments from \code{x}.

For \code{altExp<-}, a SummarizedExperiment instance with number of samples equal to \code{ncol(x)}.
Alternatively, a \code{NULL} value to remove the Experiment specified by \code{e}.
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

By default, \code{withColData=TRUE} will replace the column metadata in each alternative Experiment with that from the main object \code{x} during access to the former.
This is based on the rationale that relevant column metadata should be stored in \code{colData(x)}, given that samples are shared across the main and alternative Experiments.
Users can set this to \code{FALSE} for efficiency if they are not interested in copying the column metadata, or in the rarer case where the alternative Experiment has its own column metadata fields that are of interest.
}

\value{
\code{altExpNames} returns a character vector of names for the alternative Experiments.

\code{altExps} returns a named \linkS4class{SimpleList} of SummarizedExperiment objects.

\code{altExp} returns a SummarizedExperiment instance or a subclass thereof.

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

# Getting alternative Experimental data.
altExpNames(sce)
altExp(sce, "spike-in")
altExp(sce, 2)

# Setting alternative Experimental data.
altExpNames(sce) <- c("ERCC", "Ab")
altExp(sce, "ERCC") <- se1[1:2,]
}
