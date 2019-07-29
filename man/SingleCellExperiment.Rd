\name{SingleCellExperiment}
\alias{SingleCellExperiment}
\alias{SingleCellExperiment-class}
\alias{coerce,SummarizedExperiment,SingleCellExperiment-method}
\alias{coerce,RangedSummarizedExperiment,SingleCellExperiment-method}
\alias{updateObject,SingleCellExperiment-method}

\docType{class}

\title{SingleCellExperiment class}
\description{A description of the SingleCellExperiment class for storing single-cell sequencing data.}

\usage{
SingleCellExperiment(..., reducedDims=list(), altExps=list())
}

\arguments{
    \item{...}{Arguments to pass to the \code{\link{SummarizedExperiment}} constructor.}
    \item{reducedDims}{A list containing matrices of cell coordinates in reduced space.
    Each matrix should have the same number of rows, see \code{\link{reducedDims}}.}
    \item{altExps}{A list containing \linkS4class{SummarizedExperiment} instances of alternative feature types.
    Each element should have the same number of columns, see \code{\link{altExps}}.}
}

\value{
A SingleCellExperiment object is returned from the constructor.
}

\details{
The SingleCellExperiment class inherits from the \linkS4class{SummarizedExperiment} class.
It provides additional support for storing dimensionality reduction results (e.g., PCA, t-SNE) via \code{\link{reducedDims}}.
It also supports the storage of alternative feature types (e.g., spike-ins) via the \code{\link{altExps}}.

The intended use of this class is the same as that for SummarizedExperiment instances.
Rows should represent genomic features such as genes, while columns represent samples - in this case, single cells.
Different quantifications (e.g., counts, CPMs, log-expression) can be stored simultaneously in the \code{\link{assays}} slot.
Row and column metadata can be attached using \code{\link{rowData}} and \code{\link{colData}}, respectively.

The SingleCellExperiment contains internal metadata slots that are not intended for external use.
Please use the appropriate getter/setter functions instead, such as \code{\link{sizeFactors}}.
Package developers should refer to the suggestions in \code{?\link{int_metadata}}.
}

\seealso{
\code{\link{sizeFactors}},
\code{\link{reducedDims}}
and \code{\link{altExps}} for specific getter and setter methods.
}

\author{
Aaron Lun and Davide Risso
}

\examples{
ncells <- 100
u <- matrix(rpois(20000, 5), ncol=ncells)
v <- log2(u + 1)

pca <- matrix(runif(ncells*5), ncells)
tsne <- matrix(rnorm(ncells*2), ncells)

sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v),
    reducedDims=SimpleList(PCA=pca, tSNE=tsne))
sce

## coercion from SummarizedExperiment
se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
as(se, "SingleCellExperiment")
}

