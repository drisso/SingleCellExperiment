\name{Combining LEMs}

\alias{rbind,LinearEmbeddingMatrix-method}
\alias{cbind,LinearEmbeddingMatrix-method}

\title{LEM combining methods}

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
Note that \code{featureLoadings} and \code{factorData} will only be taken from the first element in the list;
no checks are performed to determine whether they are consistent or not across objects.

For \code{cbind}, LinearEmbeddingMatrix objects are combined columns-wise, i.e., columns in successive objects are appended to the first object.
This corresponds to adding more factors to the first object.
\code{featureLoadings} will also be combined column-wise across objects, provided that the number of features is the same across objects.
Similarly, \code{factorData} will be combined row-wise across objects.

Combining objects with and without row names will result in the removal of all row names; similarly for column names.
Duplicate row names are currently supported by duplicate column names are not, and will be de-duplicated appropriately.
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
