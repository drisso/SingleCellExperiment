\name{Miscellaneous methods}
\alias{objectVersion} 
\alias{objectVersion,SingleCellExperiment-method}
\alias{show,SingleCellExperiment-method} 
\alias{colData,SingleCellExperiment-method}
\alias{rowData,SingleCellExperiment-method}

\title{Miscellaneous methods}
\description{Various methods for the SingleCellExperiment class.}

\usage{
\S4method{show}{SingleCellExperiment}(object)

\S4method{objectVersion}{SingleCellExperiment}(x)

\S4method{colData}{SingleCellExperiment}(x, internal=FALSE)

\S4method{rowData}{SingleCellExperiment}(x, internal=FALSE)
}

\arguments{
\item{x, object}{A SingleCellExperiment object.}
\item{internal}{Whether the information contained in the internal slots should be returned.}
}

\details{
The \code{show} method will print out information about the data contained in \code{object}.
This describes the stored assays and row/column metadata, as done in \code{\link{show,SummarizedExperiment-method}}.
The names of the \code{reducedDims} slot and the names of the spike-ins (see \code{\link{spikeNames}}) are also reported.

The \code{objectVersion} method will return the version of the package with which \code{x} was constructed.
This is useful for checking if the object is up to date relative to the latest versions of the package.

It may sometimes be useful to return both the visible and the internal \code{colData} in a single DataFrame (see \code{\link{SingleCellExperiment-class}}).
This can be achieved by using \code{colData(x, internal=TRUE)}, which will return the stored \code{colData} along with the \code{int_colData} (currently the \code{sizeFactors}).
Similarly, \code{rowData(x, internal=TRUE)} will return the stored \code{rowData} along with the \code{int_rowData} (currently the columns corresponding to \code{isSpike}).
Warnings will be raised in the unlikely event of any name clashes.
}

\value{
For \code{show}, a message is printed to screen describing the data stored in \code{object}.

For \code{objectVersion}, an object of class \code{\link{package_version}} is returned.

\code{colData} and \code{rowData} return a \code{DataFrame}.
}

\seealso{
\code{\link{spikeNames}}
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
show(sce)
objectVersion(sce)
}

