\name{Combining methods}

\alias{rbind,SingleCellExperiment-method}
\alias{cbind,SingleCellExperiment-method}

\title{Combining methods}

\description{
Methods to combine SingleCellExperiment objects.
}

\usage{
\S4method{rbind}{SingleCellExperiment}(..., deparse.level=1)

\S4method{cbind}{SingleCellExperiment}(..., deparse.level=1)
}

\arguments{
\item{...}{One or more SingleCellExperiment objects.}
\item{deparse.level}{An integer scalar; see \code{?base::\link[base]{cbind}} for a description of this argument.}
}

\details{
For \code{rbind}, SingleCellExperiment objects are combined row-wise, i.e., rows in successive objects are appended to the first object.
Column metadata, experimental metadata and \code{reducedDims} coordinates will only be taken from the first element in the list.

For \code{cbind}, SingleCellExperiment objects are combined columns-wise, i.e., columns in successive objects are appended to the first object.
\code{reducedDims} coordinates will be combined row-wise to reflect the addition or more cells.
Row and experimental metadata will only be taken from the first element in the list.
}

\value{
A SingleCellExperiment object containing all rows/columns of the supplied objects.
}

\author{
Aaron Lun
}

\examples{
example(SingleCellExperiment, echo=FALSE) # using the class example
rbind(sce, sce)
cbind(sce, sce)
dim(reducedDims(sce)[[1]])
}
