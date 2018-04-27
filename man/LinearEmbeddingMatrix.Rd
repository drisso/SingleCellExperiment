\name{LinearEmbeddingMatrix}
\alias{LinearEmbeddingMatrix}
\alias{LinearEmbeddingMatrix-class}

\docType{class}

\title{LinearEmbeddingMatrix class}
\description{A description of the LinearEmbeddingMatrix class for storing low-dimensional embeddings from linear dimensionality reduction methods.}

\usage{
LinearEmbeddingMatrix(sampleFactors = matrix(nrow = 0, ncol = 0),
    featureLoadings = matrix(nrow = 0, ncol = 0), factorData = NULL,
    metadata = list())
}

\arguments{
    \item{sampleFactors}{A matrix-like object of sample embeddings, where rows are samples and columns are factors.}
    \item{featureLoadings}{A matrix-like object of feature loadings, where rows are features and columns are factors.}
    \item{factorData}{A DataFrame containing factor-level information, with one row per factor.}
    \item{metadata}{An optional list of arbitrary content describing the overall experiment.}
}

\value{
A LinearEmbeddingMatrix object is returned from the constructor.
}

\details{
The LinearEmbeddingMatrix class is a matrix-like object that supports \code{dim}, \code{dimnames} and \code{as.matrix}.
It is designed for the storage of results from linear dimensionality reduction methods like principal components analysis (PCA),
factor analysis and non-negative matrix factorization.

The \code{sampleFactors} slot is intended to store The low-dimensional representation of the samples, such as the principal coordinates from PCA.
The feature loadings contributing to each factor are stored in \code{featureLoadings}, and should have the same number of columns as \code{sampleFactors}.
The \code{factorData} stores additional factor-level information, such as the percentage of variance explained by each factor,
and should have the same number of rows as \code{sampleFactors}.

The intended use of this class is to allow PCA and other results to be stored in the \code{\link{reducedDims}} slot of a SingleCellExperiment object.
This means that feature loadings remain attached to the embedding, allowing it to be used in downstream analyses.
}

\author{
Aaron Lun, Davide Risso and Keegan Korthauer
}

\examples{
lem <- LinearEmbeddingMatrix(matrix(rnorm(1000), ncol=5),
    matrix(runif(20000), ncol=5))
lem
}

