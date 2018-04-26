\name{Getter/setter methods}
\alias{sampleFactors,LinearEmbeddingMatrix-method} 
\alias{featureLoadings,LinearEmbeddingMatrix-method} 
\alias{factorData,LinearEmbeddingMatrix-method} 
\alias{sampleFactors<-,LinearEmbeddingMatrix-method} 
\alias{featureLoadings<-,LinearEmbeddingMatrix-method} 
\alias{factorData<-,LinearEmbeddingMatrix-method} 

\title{LinearEmbeddingMatrix getters/setters}
\description{Getter/setter methods for the LinearEmbeddingMatrix class.}

\usage{
\S4method{sampleFactors}{LinearEmbeddingMatrix}(x)

\S4method{sampleFactors}{LinearEmbeddingMatrix}(x) <- value

\S4method{featureLoadings}{LinearEmbeddingMatrix}(x)

\S4method{featureLoadings}{LinearEmbeddingMatrix}(x) <- value

\S4method{factorData}{LinearEmbeddingMatrix}(x)

\S4method{factorData}{LinearEmbeddingMatrix}(x) <- value
}

\arguments{
\item{x}{A LinearEmbeddingMatrix object.}
\item{value}{An appropriate value to assign to the relevant slot.}
}

\details{
Any \code{value} to assign to \code{sampleFactors} and \code{featureLoadings} should be matrix-like objects,
while \code{factorData} should be a DataFrame - ee \linkS4class{LinearEmbeddingMatrix} for details.
}

\value{
The value of the appropriate slot for the getter methods.

For the setter methods, a LinearEmbeddingMatrix is returned with modifications to the relevant slots.
}

\seealso{
\code{\link{LinearEmbeddingMatrix}}
}

\examples{
example(LinearEmbeddingMatrix, echo=FALSE) # Using the class example

sampleFactors(lem)
sampleFactors(lem) <- sampleFactors(lem) * -1

featureLoadings(lem)
featureLoadings(lem) <- featureLoadings(lem) * -1

factorData(lem)
factorData(lem)$whee <- 1
}
