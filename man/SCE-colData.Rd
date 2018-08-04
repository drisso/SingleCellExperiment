\name{SCE colData}
\alias{colData,SingleCellExperiment-method}
\alias{rowData,SingleCellExperiment-method}

\title{SCE col/rowData methods}
\description{Methods to obtain column- or row-level metadata from the SingleCellExperiment class.}

\usage{
\S4method{colData}{SingleCellExperiment}(x, ..., internal=FALSE)

\S4method{rowData}{SingleCellExperiment}(x, ..., internal=FALSE)
}

\arguments{
\item{x}{A SingleCellExperiment object.}
\item{...}{Further arguments to be passed to the methods for \linkS4class{SummarizedExperiment} objects, such as \code{use.names}.}
\item{internal}{Whether the information contained in the internal slots should be returned.}
}

\details{
It may sometimes be useful to return both the visible and the internal \code{colData} in a single DataFrame (see \code{\link{SingleCellExperiment-class}}).
This can be achieved by using \code{colData(x, internal=TRUE)}, which will return the stored \code{colData} along with the \code{int_colData} (currently the \code{\link{sizeFactors}}).
Similarly, \code{rowData(x, internal=TRUE)} will return the stored \code{rowData} along with the \code{int_rowData} (currently the columns corresponding to \code{\link{isSpike}}).
Warnings will be raised in the event of any name clashes.
}

\value{
\code{colData} and \code{rowData} return a \code{DataFrame} with number of rows equal to \code{ncol(x)} and \code{nrow(x)}, respectively.
}

\seealso{
\code{\link{SingleCellExperiment}}
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
sizeFactors(sce) <- runif(ncol(sce))
isSpike(sce, "ERCC") <- rbinom(nrow(sce), 1, 0.2)==1
rowData(sce, internal=TRUE)
colData(sce, internal=TRUE)
}

