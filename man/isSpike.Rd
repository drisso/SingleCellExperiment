\name{Spike-in methods}

\alias{isSpike}
\alias{isSpike,SingleCellExperiment,character-method}
\alias{isSpike,SingleCellExperiment,missing-method}
\alias{isSpike,SingleCellExperiment,NULL-method}

\alias{isSpike<-}
\alias{isSpike<-,SingleCellExperiment,character-method}
\alias{isSpike<-,SingleCellExperiment,missing-method}
\alias{isSpike<-,SingleCellExperiment,NULL-method}

\alias{spikeNames}
\alias{spikeNames,SingleCellExperiment-method}

\alias{clearSpikes}
\alias{clearSpikes,SingleCellExperiment-method}

\title{Spike-in methods}
\description{Gets or sets the rows corresponding to spike-in transcripts.}

\usage{
\S4method{isSpike}{SingleCellExperiment,character}(x, type)

\S4method{isSpike}{SingleCellExperiment,missing}(x, type)

\S4method{isSpike}{SingleCellExperiment,NULL}(x, type)

\S4method{isSpike}{SingleCellExperiment,character}(x, type) <- value

\S4method{isSpike}{SingleCellExperiment,missing}(x, type) <- value

\S4method{isSpike}{SingleCellExperiment,NULL}(x, type) <- value

\S4method{clearSpikes}{SingleCellExperiment}(x)

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
This is useful not only when different spike-ins are used, but also for different mixtures of the same set (e.g., ERCC mixes 1 and 2).
The names of all available spike-in sets can be obtained using \code{spikeNames}.

To remove spike-ins for a particular set, \code{value} should be set to \code{NULL} when using \code{isSpike<-}.
To remove all spike-in information, \code{clearSpikes} should be used to obtain a new SingleCellExperiment object with no spike-ins specified.

In previous versions (<= 1.1.1), if \code{value} was \code{NULL} in \code{isSpike<-}, all existing spike-in sets would be removed.
This behaviour is now deprecated, and \code{clearSpikes} should be used instead.
Also, if \code{type} was missing or \code{NULL} for \code{isSpike<-}, the spike-in set would be automatically assigned an empty name.
This is also deprecated, and all spike-ins should be given a user-supplied name.

The \code{isSpike} getter methods will return a logical vector indicatng which rows represent spike-ins of the set specified by \code{type}.
If \code{type} is missing or \code{NULL}, the vector will instead indicate whether each row is in \emph{any} spike-in set.
If \code{type} is specified but not available, an error will be raised.
}

\value{
For \code{isSpike}, a logical vector is returned indicating whether each row is in the specified set \code{type} or any set.

For \code{isSpike<-}, a SingleCellExperiment is returned with spike-in information stored in the internal metadata fields.

For \code{spikeNames}, a character vector is returned containing the names of available spike-in sets.

For \code{clearSpikes}, a SingleCellExperiment is returned with no spike-in information.
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
