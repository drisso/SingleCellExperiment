\name{Subsetting LEMs}

\alias{[,LinearEmbeddingMatrix,ANY-method}
\alias{[,LinearEmbeddingMatrix,ANY,ANY-method}
\alias{[,LinearEmbeddingMatrix,ANY,ANY,ANY-method}
\alias{[<-,LinearEmbeddingMatrix,ANY,ANY,LinearEmbeddingMatrix-method}

\title{LEM subsetting methods}

\description{
Methods to subset LinearEmbeddingMatrix objects.
}

\usage{
\S4method{[}{LinearEmbeddingMatrix,ANY,ANY}(x, i, j, ..., drop=TRUE)

\S4method{[}{LinearEmbeddingMatrix,ANY,ANY,LinearEmbeddingMatrix}(x, i, j) <- value
}

\arguments{
\item{x}{A LinearEmbeddingMatrix object.}
\item{i, j}{A vector of logical or integer subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.}
\item{...}{Extra arguments that are ignored.}
\item{drop}{A logical scalar indicating whether the result should be coerced to the lowest possible dimension.}
\item{value}{A LinearEmbeddingMatrix object with number of rows equal to length of \code{i} (or that of \code{x}, if \code{i} is not specified).
The number of columns must be equal to the length of \code{j} (or number of columns in \code{x}, if \code{j} is not specified).
}
}

\value{
For \code{[}, a subsetted LinearEmbeddingMatrix object is returned.

For \code{[<-}, a modified LinearEmbeddingMatrix object is returned.
}

\details{
Subsetting yields a LinearEmbeddingMatrix object containing the specified rows (samples) and columns (factors).
If column subsetting is performed, values of \code{featureLoadings} and \code{factorData} will be modified to retain only the selected factors.

If \code{drop=TRUE} and the subsetting would produce dimensions of length 1, those dimensions are dropped and a vector is returned directly from \code{sampleFactors}.
This mimics the expected behaviour from a matrix-like object.
Users should set \code{drop=FALSE} to ensure that a LinearEmbeddingMatrix is returned.

For subset replacement, if neither \code{i} or \code{j} are set, \code{x} will be effectively replaced by \code{value}.
However, row and column names will \emph{not} change, consistent with replacement in ordinary matrices.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{LinearEmbeddingMatrix-class}}
}

\examples{
example(LinearEmbeddingMatrix, echo=FALSE) # using the class example

lem[1:10,]
lem[,1:5]

lem2 <- lem
lem2[1:10,] <- lem[11:20,]
}

