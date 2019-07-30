\name{Size factor methods}

\alias{sizeFactors}
\alias{sizeFactors,SingleCellExperiment-method}
\alias{sizeFactors<-,SingleCellExperiment-method}

\alias{clearSizeFactors}
\alias{clearSizeFactors,SingleCellExperiment-method}

\alias{sizeFactorNames}
\alias{sizeFactorNames,SingleCellExperiment-method}

\title{Size factors methods}
\description{Gets or sets the size factors for all cells.}

\usage{
\S4method{sizeFactors}{SingleCellExperiment}(object, type=NULL)

\S4method{sizeFactors}{SingleCellExperiment}(object, type=NULL) <- value

\S4method{clearSizeFactors}{SingleCellExperiment}(object)

\S4method{sizeFactorNames}{SingleCellExperiment}(object)
}

\arguments{
    \item{object}{A SingleCellExperiment object.}
    \item{type}{A string specifying the \emph{type} of size factor to get or set.}
    \item{value}{A numeric vector of size factors for all cells.}
}

\details{
A size factor is a scaling factor used to divide the raw counts of a particular cell to obtain normalized expression values.
The \code{sizeFactors} methods can be used to get or set size factors for all cells.

The \code{type} argument is deprecated, as are the \code{sizeFactorNames} and \code{clearSizeFactors} functions.
}

\value{
For \code{sizeFactors}, a numeric vector is returned containing size factors for all cells.

For \code{sizeFactors<-}, a SingleCellExperiment is returned with size factors stored in the internal metadata fields.

For \code{clearSizeFactors}, a SingleCellExperiment is returned with no size factor information.

For \code{sizeFactorNames}, a character vector is returned containing the names of all named size factor sets.
}

\seealso{
\code{\link{SingleCellExperiment-class}}
}

\author{
Aaron Lun
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
sizeFactors(sce) <- runif(ncol(sce))
sizeFactors(sce)
}
