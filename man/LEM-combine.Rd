\name{Combining methods}

\alias{rbind,LinearEmbeddingMatrix-method}
\alias{cbind,LinearEmbeddingMatrix-method}

\title{Combining methods}

\description{
Methods to combine LinearEmbeddingMatrix objects.
}

\usage{
\S4method{rbind}{LinearEmbeddingMatrix}(..., deparse.level=1)

\S4method{cbind}{LinearEmbeddingMatrix}(..., deparse.level=1)
}

\arguments{
\item{...}{One or more LinearEmbeddingMatrix objects.}
\item{deparse.level}{An integer scalar; see \code{?base::\link[base]{cbind}} for a description of this argument.}
}

\details{
For \code{rbind}, LinearEmbeddingMatrix objects are combined row-wise, i.e., rows in successive objects are appended to the first object.
This corresponds to adding more samples to the first object.
\code{featureLoadings} and \code{factorData} will only be taken from the first element in the list.
Matrix-like objects that are \code{not} LinearEmbeddingMatrix objects can also be used in \code{...}, 
    provided that they are not the first element and that they support \code{rbind}.

For \code{cbind}, LinearEmbeddingMatrix objects are combined columns-wise, i.e., columns in successive objects are appended to the first object.
This corresponds to adding more factors to the first object.
\code{featureLoadings} will also be combined column-wise across objects, provided that the number of features is the same across objects.
Similarly, \code{factorData} will be combined row-wise across objects.
}

\value{
A LinearEmbeddingMatrix object containing all rows/columns of the supplied objects.
}

\author{
Aaron Lun
}

\examples{
example(LinearEmbeddingMatrix, echo=FALSE) # using the class example
rbind(lem, lem)
cbind(lem, lem)
}
