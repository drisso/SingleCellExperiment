\name{SingleCellExperiment}
\alias{SingleCellExperiment}
\alias{SingleCellExperiment-class}
\alias{coerce,SummarizedExperiment,SingleCellExperiment-method}
\alias{colData,SingleCellExperiment-method}
\alias{rowData,SingleCellExperiment-method}

\docType{class}

\title{SingleCellExperiment class}
\description{A description of the SingleCellExperiment class for storing single-cell sequencing data.}

\usage{
SingleCellExperiment(..., reducedDims=SimpleList())

\S4method{colData}{SingleCellExperiment}(x, internal=FALSE)
\S4method{rowData}{SingleCellExperiment}(x, internal=FALSE)
}

\arguments{
    \item{...}{Arguments to pass to the \code{SummarizedExperiment} constructor.}
    \item{reducedDims}{A SimpleList object containing matrices of cell coordinates in reduced space.}
    \item{x}{A \code{SingleCellExperiment} object.}
    \item{internal}{Whether the information contained in the internal slots should be returned. See details.}
}

\value{
A SingleCellExperiment object is returned from the constructor.

\code{colData} and \code{rowData} return a \code{DataFrame}.
}

\details{
The SingleCellExperiment class inherits from the \code{\link{SummarizedExperiment}} class, with several additional slots:
\describe{
\item{\code{reducedDims}:}{A SimpleList containing matrices of cell coordinates.}
\item{\code{int_elementMetadata}:}{A DataFrame containing internal row metadata (for each genomic feature).}
\item{\code{int_colData}:}{A DataFrame containing internal column metadata (for each cell).}
\item{\code{int_metadata}:}{A list containing internal experiment metadata.}
}

The intended use of this class is the same as that for SummarizedExperiment instances.
Rows should represent genomic features such as genes, while columns represent samples - in this case, single cells.
Different quantifications (e.g., counts, CPMs, log-expression) can be stored simultaneously in the \code{\link{assays}} slot.
Row and column metadata can be attached using \code{\link{rowData}} and \code{\link{colData}}, respectively.

The additional \code{reducedDims} slot allows storage of results from multiple dimensionality reduction methods, e.g., PCA or t-SNE.
Each element of the SimpleList should be a matrix of coordinates for all cells from one reduction method.
The number of rows of each matrix should be equal to the number of cells in the SingleCellExperiment object.

The internal metadata slots are not intended for external use.
Please use the appropriate getter/setter functions instead, such as \code{isSpike} or \code{sizeFactors}.

It may sometimes be useful to return both the visible and the internal \code{colData} in a single \code{DataFrame}. To do so, please use \code{colData(x, internal=TRUE)}: this will return the stored \code{colData}, along with the \code{int_colData} (currently the \code{sizeFactors}).

Analogously, \code{rowData(x, internal=TRUE)} will return the stored \code{rowData} along with the \code{int_rowData} (currently the columns corresponding to \code{isSpike}).
}

\seealso{
\code{\link{isSpike}},
\code{\link{sizeFactors}},
\code{\link{reducedDims}}
}

\author{
Aaron Lun and Davide Risso
}

\examples{
ncells <- 100
v <- matrix(rnorm(20000), ncol=ncells)
u <- matrix(rpois(20000, 5), ncol=ncells)

pca <- matrix(runif(ncells*5), ncells)
tsne <- matrix(rnorm(ncells*2), ncells)

sce <- SingleCellExperiment(assays=list(counts=u, exprs=v),
    reducedDims=SimpleList(PCA=pca, tSNE=tsne))
sce

## coercion from SummarizedExperiment
se <- SummarizedExperiment(assays=list(counts=u, exprs=v))
as(se, "SingleCellExperiment")
}

