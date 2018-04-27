\name{Miscellaneous LEM}
\alias{show,LinearEmbeddingMatrix-method} 

\title{Miscellaneous LEM methods}
\description{Various methods for the LinearEmbeddingMatrix class.}

\usage{
\S4method{show}{LinearEmbeddingMatrix}(object)
}

\arguments{
\item{object}{A LinearEmbeddingMatrix object.}
}

\details{
The \code{show} method will print out information about the data contained in \code{object}.
This includes the number of samples, the number of factors, the number of genes and the fields available in \code{factorData}.
}

\value{
A message is printed to screen describing the data stored in \code{object}.
}

\seealso{
\code{\link{LinearEmbeddingMatrix}}
}

\author{
Davide Risso
}

\examples{
example(LinearEmbeddingMatrix, echo=FALSE) # Using the class example
show(lem)
}

