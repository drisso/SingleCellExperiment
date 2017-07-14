\name{Subsetting methods}

\alias{[,SingleCellExperiment,ANY-method}
\alias{[,SingleCellExperiment,ANY,ANY-method}
\alias{[,SingleCellExperiment,ANY,ANY,ANY-method}
\alias{subset,SingleCellExperiment-method}
\alias{[<-,SingleCellExperiment,ANY,ANY,SingleCellExperiment-method}

\title{Subsetting methods}

\description{
Methods to subset SingleCellExperiment objects.
}

\usage{
\S4method{[}{SingleCellExperiment,ANY,ANY}(x, i, j, ..., drop=TRUE)

\S4method{subset}{SingleCellExperiment}(x, i, j)

\S4method{[}{SingleCellExperiment,ANY,ANY,SingleCellExperiment}(x, i, j) <- value
}

\arguments{
\item{x}{A SingleCellExperiment object.}
\item{i, j}{A vector of logical or integer subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.}
\item{...}{Extra arguments to be passed to \code{\link{[,SummarizedExperiment-method}}.}
\item{drop}{A logical scalar that is ignored.}
\item{value}{A SingleCellExperiment object with number of rows equal to length of \code{i} (or that of \code{x}, if \code{i} is not specified).
The number of columns must be equal to the length of \code{j} (or number of columns in \code{x}, if \code{j} is not specified).
}
}

\value{
For \code{[} and \code{subset}, a subsetted SingleCellExperiment object is returned.

For \code{[<-}, a modified SingleCellExperiment object is returned.
}

\details{
Subsetting yields a SingleCellExperiment object containing the specified rows (features) and columns (cells).
Internal row and column metadata fields will also be subsetted so that methods such as \code{\link{isSpike}} are still valid.
If column subsetting is performed, values of the \code{reducedDims} will be modified to retain only the selected cells.

Subset assignment will replace the assay values and metadata of the specified rows or columns in \code{x} with those in \code{value}.
If both \code{i} and \code{j} are set, the relevant block of assay values will be replaced, along with the metadata for the affected rows and columns.
If neither \code{i} or \code{j} are set, \code{x} will be turned into \code{value}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{SingleCellExperiment-class}}
}

\examples{
example(SingleCellExperiment, echo=FALSE) # using the class example

sce[1:10,]
sce[,1:5]

sce2 <- sce
sce2[1:10,] <- sce[11:20,]

# Can also use subset()
subset(sce, 1, 1)

# Can also use split()
split(sce, sample(LETTERS, nrow(sce), replace=TRUE))
}

