\name{Spike-in methods}

\alias{isSpike}
\alias{isSpike,SingleCellExperiment,character-method}
\alias{isSpike,SingleCellExperiment,missing-method}
\alias{isSpike,SingleCellExperiment,NULL-method}

\alias{isSpike<-}
\alias{isSpike<-,SingleCellExperiment,character-method}

\alias{spikeNames}
\alias{spikeNames,SingleCellExperiment-method}

\title{Spike-in methods}
\description{Gets or sets the rows corresponding to spike-in transcripts.}

\usage{
\S4method{isSpike}{SingleCellExperiment,character}(x, type)

\S4method{isSpike}{SingleCellExperiment,missing}(x, type)

\S4method{isSpike}{SingleCellExperiment,NULL}(x, type)

\S4method{isSpike}{SingleCellExperiment,character}(x, type) <- value

\S4method{spikeNames}{SingleCellExperiment}(x)
}

\arguments{
    \item{x}{A SingleCellExperiment object.}
    \item{type}{A string containing the name of the spike-in set.}
    \item{value}{A vector indicating which rows correspond to spike-in transcripts.}
}

\details{
Spike-in transcripts may be added during library preparation in single-cell RNA sequencing experiments.
These usually need to be handled differently during data analysis, compared to the endogenous genes.
Thus, it is important to indicate which rows correspond to spike-in transcripts.

The \code{isSpike<-} method accepts any \code{value} that indicates which rows correspond to spike-ins.
This can be a logical or integer subsetting vector, or a vector of row names.
The \code{type} should be set to the name of the spike-in set, e.g., \code{"ERCC"} or \code{"SIRV"}.

In this manner, multiple types of spike-in sets are supported for a single experiment.
The names of all available spike-in sets can be obtained using \code{spikeNames}.
To remove spike-ins for a particular set, \code{value} should be set to \code{NULL} when using \code{isSpike<-}.

The \code{isSpike} getter methods will return a logical vector indicatng which rows represent spike-ins of the set specified by \code{type}.
If \code{type} is missing or \code{NULL}, the vector will instead indicate whether each row is in \emph{any} spike-in set.
If \code{type} is specified but not available, an error will be raised.
}

\value{
For \code{isSpike}, a logical vector is returned indicating whether each row is in the specified set \code{type} or any set.

For \code{isSpike<-}, a SingleCellExperiment is returned with spike-in information stored in the internal metadata fields.

For \code{spikeNames}, a character vector is returned containing the names of available spike-in sets.
}

\seealso{
\code{\link{SingleCellExperiment-class}}
}

\author{
Aaron Lun
}

\examples{
example(SingleCellExperiment, echo=FALSE) # Using the class example
isSpike(sce, "ERCC") <- 1:10
isSpike(sce)

isSpike(sce, "SIRV") <- 11:20
spikeNames(sce)
which(isSpike(sce))
which(isSpike(sce, "SIRV"))

isSpike(sce, "ERCC") <- NULL
spikeNames(sce)
}
