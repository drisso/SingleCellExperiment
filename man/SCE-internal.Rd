\name{SCE internals}

\alias{int_colData}
\alias{int_elementMetadata}
\alias{int_metadata}

\alias{int_colData,SingleCellExperiment-method}
\alias{int_elementMetadata,SingleCellExperiment-method}
\alias{int_metadata,SingleCellExperiment-method}

\alias{int_colData<-}
\alias{int_elementMetadata<-}
\alias{int_metadata<-}

\alias{int_colData<-,SingleCellExperiment-method}
\alias{int_elementMetadata<-,SingleCellExperiment-method}
\alias{int_metadata<-,SingleCellExperiment-method}

\title{SCE internal methods}
\description{Methods to get or set internal fields from the SingleCellExperiment class.}

\usage{
\S4method{int_elementMetadata}{SingleCellExperiment}(x)

\S4method{int_elementMetadata}{SingleCellExperiment}(x) <- value

\S4method{int_colData}{SingleCellExperiment}(x)

\S4method{int_colData}{SingleCellExperiment}(x) <- value

\S4method{int_metadata}{SingleCellExperiment}(x)

\S4method{int_metadata}{SingleCellExperiment}(x) <- value
}

\arguments{
\item{x}{A SingleCellExperiment object.}
\item{value}{For \code{int_elementMetadata}, a \linkS4class{DataFrame} with number of rows equal to \code{nrow(x)}.

For \code{int_colData}, a DataFrame with number of rows equal to \code{ncol(x)}.

For \code{int_metadata}, a list.}
}

\details{
These functions are intended for package developers who want to add protected fields to a SingleCellExperiment.
They should \emph{not} be used by ordinary users of the \pkg{SingleCellExperiment} package.

Package developers intending to use these methods should read the development vignette for instructions. 
}

\value{
A SingleCellExperiment object equivalent to \code{x} but with modified internal fields.
}

\seealso{
\code{\link{SingleCellExperiment}}
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
int_metadata(sce)$whee <- 1
}

