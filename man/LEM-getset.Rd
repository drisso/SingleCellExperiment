\name{Getter/setter methods}
\alias{sampleFactors} 
\alias{featureLoadings} 
\alias{factorData} 

\alias{sampleFactors,LinearEmbeddingMatrix-method} 
\alias{featureLoadings,LinearEmbeddingMatrix-method} 
\alias{factorData,LinearEmbeddingMatrix-method} 

\alias{sampleFactors<-} 
\alias{featureLoadings<-} 
\alias{factorData<-} 

\alias{sampleFactors<-,LinearEmbeddingMatrix-method} 
\alias{featureLoadings<-,LinearEmbeddingMatrix-method} 
\alias{factorData<-,LinearEmbeddingMatrix-method} 

\alias{as.matrix,LinearEmbeddingMatrix-method} 
\alias{dim,LinearEmbeddingMatrix-method} 
\alias{dimnames,LinearEmbeddingMatrix-method} 
\alias{dimnames<-,LinearEmbeddingMatrix-method} 
\alias{dimnames<-,LinearEmbeddingMatrix,ANY-method} 

\alias{$,LinearEmbeddingMatrix-method}
\alias{$<-,LinearEmbeddingMatrix-method}

\title{LinearEmbeddingMatrix getters/setters}
\description{Getter/setter methods for the LinearEmbeddingMatrix class.}

\usage{
\S4method{sampleFactors}{LinearEmbeddingMatrix}(x, withDimnames=TRUE)

\S4method{sampleFactors}{LinearEmbeddingMatrix}(x) <- value

\S4method{featureLoadings}{LinearEmbeddingMatrix}(x, withDimnames=TRUE)

\S4method{featureLoadings}{LinearEmbeddingMatrix}(x) <- value

\S4method{factorData}{LinearEmbeddingMatrix}(x)

\S4method{factorData}{LinearEmbeddingMatrix}(x) <- value

\S4method{as.matrix}{LinearEmbeddingMatrix}(x, ...)

\S4method{dim}{LinearEmbeddingMatrix}(x)

\S4method{dimnames}{LinearEmbeddingMatrix}(x)

\S4method{dimnames}{LinearEmbeddingMatrix}(x) <- value

\S4method{$}{LinearEmbeddingMatrix}(x, name)

\S4method{$}{LinearEmbeddingMatrix}(x, name) <- value
}

\arguments{
\item{x}{A LinearEmbeddingMatrix object.}
\item{value}{An appropriate value to assign to the relevant slot.}
\item{withDimnames}{A logical scalar indicating whether dimension names should be attached to the returned object.}
\item{name}{A string specifying a field of the \code{factorData} slot.}
\item{...}{Further arguments, ignored.}
}

\details{
Any \code{value} to assign to \code{sampleFactors} and \code{featureLoadings} should be matrix-like objects,
while \code{factorData} should be a DataFrame - ee \linkS4class{LinearEmbeddingMatrix} for details.

The \code{as.matrix} method will return the matrix of sample factors, consistent with the fact that the LinearEmbeddingMatrix mimics a sample-factor matrix.
However, unlike the \code{sampleFactors} method, this is always guaranteed to return an ordinary R matrix, even if an alternative representation was stored in the slot.
This ensures consistency with \code{as.matrix} methods for other matrix-like S4 classes. 

For assignment to \code{dimnames}, a list of length 2 should be used containing vectors of row and column names.
}

\value{
For the getter methods \code{sampleFactors}, \code{featureLoadings} and \code{factorData}, the value of the slot with the same name is returned.
For the corresponding setter methods, a LinearEmbeddingMatrix is returned with modifications to the named slot.

For \code{dim}, the dimensions of the \code{sampleFactors} slot are returned in an integer vector of length 2.
For \code{dimnames}, a list of length 2 containing the row and column names is returned.
For \code{as.matrix}, an ordinary matrix derived from \code{sampleFactors} is returned.

For \code{$}, the value of the named field of the \code{factorData} slot is returned.
For \code{$<-}, a LinearEmbeddingMatrix is returned with the modified field in \code{factorData}.
}

\seealso{
\code{\link{LinearEmbeddingMatrix}}
}

\author{
Keegan Korthauer, Davide Risso and Aaron Lun
}

\examples{
example(LinearEmbeddingMatrix, echo=FALSE) # Using the class example

sampleFactors(lem)
sampleFactors(lem) <- sampleFactors(lem) * -1

featureLoadings(lem)
featureLoadings(lem) <- featureLoadings(lem) * -1

factorData(lem)
factorData(lem)$whee <- 1

nrow(lem)
ncol(lem)
colnames(lem) <- LETTERS[seq_len(ncol(lem))]
as.matrix(lem)
}
