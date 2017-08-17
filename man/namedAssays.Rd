\name{namedAssays}
\alias{counts}
\alias{counts<-}
\alias{counts,SingleCellExperiment-method}
\alias{counts<-,SingleCellExperiment-method}

\alias{normcounts}
\alias{normcounts<-}
\alias{normcounts,SingleCellExperiment-method}
\alias{normcounts<-,SingleCellExperiment-method}

\alias{logcounts}
\alias{logcounts<-}
\alias{logcounts,SingleCellExperiment-method}
\alias{logcounts<-,SingleCellExperiment-method}

\alias{cpm}
\alias{cpm<-}
\alias{cpm,SingleCellExperiment-method}
\alias{cpm<-,SingleCellExperiment-method}

\alias{tpm}
\alias{tpm<-}
\alias{tpm,SingleCellExperiment-method}
\alias{tpm<-,SingleCellExperiment-method}

\title{Named assay fields}
\description{Convenience methods to get or set named assay fields.}

\usage{
\S4method{counts}{SingleCellExperiment}(object)
\S4method{counts}{SingleCellExperiment}(object) <- value

\S4method{normcounts}{SingleCellExperiment}(object)
\S4method{normcounts}{SingleCellExperiment}(object) <- value

\S4method{logcounts}{SingleCellExperiment}(object)
\S4method{logcounts}{SingleCellExperiment}(object) <- value

\S4method{cpm}{SingleCellExperiment}(object)
\S4method{cpm}{SingleCellExperiment}(object) <- value

\S4method{tpm}{SingleCellExperiment}(object)
\S4method{tpm}{SingleCellExperiment}(object) <- value
}

\arguments{
\item{object}{A SingleCellExperiment object.}
\item{value}{A numeric matrix of the same dimensions as \code{object}.}
}

\details{
These are wrapper methods for getting or setting \code{assay(object, i=X)} where \code{X} is the name of the method.
For example, \code{counts} will get or set \code{X="counts"}.
This provide some convenience for users as well as encouraging standardization of naming across packages.

Our suggested interpretation of the fields are as follows:
\describe{
\item{\code{counts}:}{Raw count data, e.g., number of reads or transcripts.}
\item{\code{normcounts}:}{Normalized values on the same scale as the original counts.
For example, counts divided by cell-specific size factors that are centred at unity.}
\item{\code{logcounts}:}{Log-transformed counts or count-like values.
In most cases, this will be defined as log-transformed \code{normcounts}, e.g., using log base 2 and a pseudo-count of 1.}
\item{\code{cpm}:}{Counts-per-million.
This is the read count for each gene in each cell, divided by the library size of each cell in millions.}
\item{\code{tpm}:}{Transcripts-per-million.
This is the number of transcripts for each gene in each cell, divided by the total number of transcripts in that cell (in millions).}
}
}

\value{
Each method returns a matrix from the correspondingly named field in the \code{assays} slot.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{SingleCellExperiment}}
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
counts(sce) <- matrix(rnorm(nrow(sce)*ncol(sce)), ncol=ncol(sce))
dim(counts(sce))

# One possible way of computing normalized "counts"
sf <- 2^rnorm(ncol(sce))
sf <- sf/mean(sf) 
normcounts(sce) <- t(t(counts(sce))/sf)
dim(normcounts(sce))

# One possible way of computing log-counts
logcounts(sce) <- log2(normcounts(sce)+1)
dim(normcounts(sce))
}
