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
\alias{reducedDimNames<-}
\alias{reducedDimNames<-,SingleCellExperiment,character-method}

\title{Reduced dimensions methods}
\description{Methods to get or set the dimensionality reduction results.}

\usage{
\S4method{reducedDim}{SingleCellExperiment}(x, type=1, withDimnames=TRUE)

\S4method{reducedDim}{SingleCellExperiment}(x, type=1, ...) <- value

\S4method{reducedDims}{SingleCellExperiment}(x, withDimnames=TRUE)

\S4method{reducedDims}{SingleCellExperiment}(x) <- value

\S4method{reducedDimNames}{SingleCellExperiment}(x)

\S4method{reducedDimNames}{SingleCellExperiment,character}(x) <- value
}

\arguments{
\item{x}{A \linkS4class{SingleCellExperiment} object.}
\item{type}{A string containing the name for the dimensionality reduction results or a numeric index containing the position of the desired dimenionality reduction result.}
\item{withDimnames}{A logical scalar indicating whether each set of results should be returned with row names matching \code{colnames(x)}.
This can be turned off for greater efficiency.}
\item{...}{Additional arguments that are ignored.}
\item{value}{
For \code{reducedDim<-}, a matrix (usually double-precision) of coordinates, for each cell (row) and dimension (column).

For \code{reducedDims<-}, a named SimpleList object containing such matrices.
Each entry should have the same number of rows corresponding to the set of cells in \code{x}.

For \code{reducedDimNames<-}, a character vector containing the names of all dimensionality reduction results.}
}

\details{
Dimensionality reduction is often used to interpreting the results of single-cell data analysis.
These methods allow the results of dimensionality reduction methods to be stored in a SingleCellExperiment object.
Multiple results can be stored in a single object by assigning to different \code{type} in \code{reducedDim<-}.

If \code{value} is \code{NULL} for \code{reducedDim<-}, the set of results corresponding to \code{type} is removed from the object.
If \code{value} is \code{NULL} for \code{reducedDims<-}, all dimensionality reduction results are removed.

Unnamed results assigned via \code{reducedDim<-} or \code{reducedDims<-} will be assigned empty names.
}

\value{
For \code{reducedDim}, a numeric matrix containing coordinates for cells (rows) and dimensions (columns).

For \code{reducedDims}, a named \linkS4class{SimpleList} of matrices containing one or more dimensionality reduction results.

For \code{reducedDimNames}, a character vector containing the names of the elements in \code{reducedDims}.

For the replacement methods, a SingleCellExperiment object with updated dimensionality reduction results.
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
