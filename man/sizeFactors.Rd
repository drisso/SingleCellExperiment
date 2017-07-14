\name{Size factor methods}
\alias{sizeFactors}
\alias{sizeFactors,SingleCellExperiment-method}
\alias{sizeFactors<-,SingleCellExperiment-method}

\title{Size factors methods}
\description{Gets or sets the size factors for all cells.}

\usage{
\S4method{sizeFactors}{SingleCellExperiment}(object, type=NULL)

\S4method{sizeFactors}{SingleCellExperiment}(object, type=NULL) <- value
}

\arguments{
    \item{object}{A SingleCellExperiment object.}
    \item{type}{A string specifying the \emph{type} of size factor to get or set.}
    \item{value}{A numeric vector of size factors for all cells.}
}

\details{
A size factor is a scaling factor used to divide the raw counts of a particular cell to obtain normalized expression values.
The \code{sizeFactors} methods can be used to get or set size factors for all cells.

The \code{type} argument allows storage of multiple vectors of size factors (e.g., different values for spike-ins versus endogenous genes).
If \code{type} is \code{NULL}, a \dQuote{default} set of size factors is stored or returned.

If \code{value} is \code{NULL} for \code{isSpike<-}, size factors of \code{type} will be removed from \code{object}.
}

\value{
For \code{sizeFactors}, a numeric vector is returned containing size factors of the set \code{type} for all cells.
If \code{type} is not available, \code{NULL} is returned instead.

For \code{sizeFactors<-}, a SingleCellExperiment is returned with size factors stored in the internal metadata fields.
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

sizeFactors(sce, "ERCC") <- runif(ncol(sce)) 
sizeFactors(sce, "ERCC")
sizeFactors(sce) # unchanged.

sizeFactors(sce, "ERCC") <- NULL
sizeFactors(sce, "ERCC") 
}
