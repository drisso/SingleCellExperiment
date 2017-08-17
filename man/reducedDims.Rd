\name{Reduced dimensions}

\alias{reducedDim}
\alias{reducedDim,SingleCellExperiment-method}

\alias{reducedDims}
\alias{reducedDims,SingleCellExperiment-method}

\alias{reducedDim<-}
\alias{reducedDim<-,SingleCellExperiment-method}
\alias{reducedDims<-}
\alias{reducedDims<-,SingleCellExperiment-method}

\alias{reducedDimNames}
\alias{reducedDimNames,SingleCellExperiment-method}

\title{Reduced dimensions methods}
\description{Methods to get or set the dimensionality reduction results.}

\usage{
\S4method{reducedDim}{SingleCellExperiment}(x, type)

\S4method{reducedDim}{SingleCellExperiment}(x, type) <- value

\S4method{reducedDims}{SingleCellExperiment}(x)

\S4method{reducedDims}{SingleCellExperiment}(x) <- value
}

\arguments{
\item{x}{A SingleCellExperiment object.}
\item{type}{A string containing the name for the dimensionality reduction results or a numeric index containing the position of the desired dimenionality reduction result.}
\item{value}{
For \code{reducedDim<-}, a matrix (usually double-precision) of coordinates, for each cell (row) and dimension (column).
For \code{reducedDims<-}, a named SimpleList object containing such matrices.}
}

\details{
Dimensionality reduction is often used to interpreting the results of single-cell data analysis.
These methods allow the results of dimensionality reduction methods to be stored in a SingleCellExperiment object.
Multiple results can be stored in a single object by assigning to different \code{type} in \code{reducedDim<-}.

If \code{type} is \code{NULL} or missing for \code{reducedDim}, the first set of dimensionality reduction results is returned (or \code{NULL}, if no results are present).
If \code{value} is \code{NULL} for \code{reducedDim<-}, the set of results corresponding to \code{type} is removed from the object.
}

\value{
For \code{reducedDim}, a numeric matrix is returned containing coordinates for cells (rows) and dimensions (columns).

For \code{reducedDims}, a named SimpleList of matrices is returned, with one matrix for each type of dimensionality reduction method.

For \code{reducedDim<-} and \code{reducedDims<-}, a SingleCellExperiment object is returned with updated results in the \code{reducedDims} slot.

For \code{reducedDimNames}, a character vector containing the names of the elements in \code{reducedDims}.
}

\author{
Aaron Lun
}

\seealso{
    \code{\link{SingleCellExperiment-class}}
}

\examples{
example(SingleCellExperiment, echo=FALSE)
reducedDim(sce, "PCA")
reducedDim(sce, "tSNE")
reducedDims(sce)

reducedDim(sce, "PCA") <- NULL
reducedDims(sce)

reducedDims(sce) <- SimpleList()
reducedDims(sce)
}
