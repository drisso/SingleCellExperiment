#' Alternative Experiment methods
#'
#' @description
#' In some experiments, different features must be normalized differently or have different row-level metadata.
#' Typical examples would be for spike-in transcripts in plate-based experiments and antibody or CRISPR tags in CITE-seq experiments.
#' These data cannot be stored in the main \code{assays} of the \linkS4class{SingleCellExperiment} itself.
#' However, it is still desirable to store these features \emph{somewhere} in the SingleCellExperiment.
#' This simplifies book-keeping in long workflows and ensure that samples remain synchronised.
#'
#' To facilitate this, the \linkS4class{SingleCellExperiment} class allows for \dQuote{alternative Experiments}.
#' Nested \linkS4class{SummarizedExperiment}-class objects are stored inside the SingleCellExperiment object \code{x}, in a manner that guarantees that the nested objects have the same columns in the same order as those in \code{x}.
#' Methods are provided to enable convenient access to and manipulation of these alternative Experiments.
#' Each alternative Experiment should contain experimental data and row metadata for a distinct set of features.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{altExp(x, e, withDimnames=TRUE, withColData=FALSE)}:}{
#' Retrieves a \linkS4class{SummarizedExperiment} containing alternative features (rows) for all cells (columns) in \code{x}.
#' \code{e} should either be a string specifying the name of the alternative Experiment in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired Experiment, defaulting to the first Experiment is missing.
#'
#' If \code{withDimnames=TRUE}, the column names of the output object are set to \code{colnames(x)}.
#' In addition, if \code{withColData=TRUE}, \code{\link{colData}(x)} is \code{cbind}ed to the front of the column data of the output object.
#' }
#' \item{\code{altExpNames(x)}:}{
#' Returns a character vector containing the names of all alternative Experiments in \code{x}.
#' This is guaranteed to be of the same length as the number of results, though the names may not be unique.
#' }
#' \item{\code{altExps(x, withDimnames=TRUE, withColData=FALSE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more \linkS4class{SummarizedExperiment} objects.
#' Each object is guaranteed to have the same number of columns, in a 1:1 correspondence to those in \code{x}.
#'
#' If \code{withDimnames=TRUE}, the column names of each output object are set to \code{colnames(x)}.
#' In addition, if \code{withColData=TRUE}, \code{\link{colData}(x)} is \code{cbind}ed to the front of the column data of each output object.
#' }
#' }
#'
#' @section Single-object setter:
#' \code{altExp(x, e, withDimnames=TRUE, withColData=FALSE) <- value} will add or replace an alternative Experiment
#' in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{e} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{e} is missing, \code{value} is assigned to the first result.
#' If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed1"}.
#' \item If \code{e} is a numeric scalar, it must be within the range of existing results, and \code{value} will be assigned to the result at that index.
#' \item If \code{e} is a string and a result exists with this name, \code{value} is assigned to to that result.
#' Otherwise a new result with this name is append to the existing list of results.
#' }
#'
#' \code{value} is expected to be a SummarizedExperiment object with number of columns equal to \code{ncol(x)}.
#' Alternatively, if \code{value} is \code{NULL}, the alternative Experiment at \code{e} is removed from the object.
#'
#' If \code{withDimnames=TRUE}, the column names of \code{value} are checked against those of \code{x}.
#' A warning is raised if these are not identical, with the only exception being when \code{value=NULL}.
#' This is inspired by the argument of the same name in \code{\link{assay<-}}.
#'
#' If \code{withColData=TRUE}, we assume that the left-most columns of \code{colData(value)} are identical to \code{colData(x)}.
#' If so, these columns are removed, effectively reversing the \code{withColData=TRUE} setting for the \code{altExp} getter.
#' Otherwise, a warning is raised.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{altExps(x, withDimnames=TRUE, withColData=FALSE) <- value}:}{
#' Replaces all alterrnative Experiments in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of SummarizedExperiment objects
#' with number of columns equal to \code{ncol(x)}.
#'
#' If \code{value} is named, those names will be used to name the alternative Experiments in \code{x}.
#' Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed"}.
#'
#' If \code{value} is \code{NULL}, all alternative Experiments in \code{x} are removed.
#'
#' If \code{value} is a \linkS4class{Annotated} object, any \code{\link{metadata}} will be retained in \code{altExps(x)}.
#' If \code{value} is a \linkS4class{Vector} object, any \code{\link{mcols}} will also be retained.
#'
#' If \code{withDimnames=TRUE}, the column names of each entry of \code{value} are checked against those of \code{x}.
#' A warning is raised if these are not identical.
#'
#' If \code{withColData=TRUE}, we assume that the left-most columns of the \code{colData} for each entry of \code{value} are identical to \code{colData(x)}.
#' If so, these columns are removed, effectively reversing the \code{withColData=TRUE} setting for the \code{altExps} getter.
#' Otherwise, a warning is raised.
#' }
#' \item{\code{altExpNames(x) <- value}:}{
#' Replaces all names for alternative Experiments in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of results currently in \code{x}.
#' }
#' }
#'
#' \code{removeAltExps(x)} will remove all alternative Experiments from \code{x}.
#' This has the same effect as \code{altExps(x) <- NULL} but may be more convenient as it directly returns a SingleCellExperiment.
#'
#' @section Main Experiment naming:
#' The alternative Experiments are naturally associated with names (\code{e} during assignment).
#' However, we can also name the main Experiment in a \linkS4class{SingleCellExperiment} \code{x}:
#' \describe{
#' \item{\code{mainExpName(x) <- value}:}{
#' Set the name of the main Experiment to a non-\code{NA} string \code{value}.
#' This can also be used to unset the name if \code{value=NULL}.
#' }
#' \item{\code{mainExpName(x)}:}{
#' Returns a string containing the name of the main Experiment.
#' This may also be \code{NULL} if no name is specified.
#' }
#' }
#' The presence of a non-\code{NULL} main Experiment name is helpful for functions like \code{\link{swapAltExp}}.
#' An appropriate name is automatically added by functions like \code{\link{splitAltExps}}.
#' 
#' Note that, if a SingleCellExperiment is assigned as an alternative Experiment to another SingleCellExperiment via \code{altExp(x, e) <- value},
#' no attempt is made to synchronize \code{mainExpName(value)} with \code{e}.
#' In such cases, we suggest setting \code{mainExpName(value)} to \code{NULL} to avoid any confusion during interpretation.
#' 
#' @seealso
#' \code{\link{splitAltExps}}, for a convenient way of adding alternative Experiments from existing features.
#'
#' \code{\link{swapAltExp}}, to swap the main and alternative Experiments.
#'
#' @author Aaron Lun
#'
#' @examples
#' example(SingleCellExperiment, echo=FALSE) # Using the class example
#' dim(counts(sce))
#'
#' # Mocking up some alternative Experiments.
#' se1 <- SummarizedExperiment(matrix(rpois(1000, 5), ncol=ncol(se)))
#' rowData(se1)$stuff <- sample(LETTERS, nrow(se1), replace=TRUE)
#' se2 <- SummarizedExperiment(matrix(rpois(500, 5), ncol=ncol(se)))
#' rowData(se2)$blah <- sample(letters, nrow(se2), replace=TRUE)
#'
#' # Setting the alternative Experiments.
#' altExp(sce, "spike-in") <- se1
#' altExp(sce, "CRISPR") <- se2
#'
#' # Getting alternative Experimental data.
#' altExpNames(sce)
#' altExp(sce, "spike-in")
#' altExp(sce, 2)
#'
#' # Setting alternative Experimental data.
#' altExpNames(sce) <- c("ERCC", "Ab")
#' altExp(sce, "ERCC") <- se1[1:2,]
#'
#' @name altExps
#' @docType methods
#' @aliases
#' altExp altExps altExpNames
#' altExp,SingleCellExperiment,missing-method
#' altExp,SingleCellExperiment,numeric-method
#' altExp,SingleCellExperiment,character-method
#' altExps,SingleCellExperiment-method
#' altExpNames,SingleCellExperiment-method
#' altExp<- altExps<- altExpNames<-
#' altExp<-,SingleCellExperiment,missing-method
#' altExp<-,SingleCellExperiment,numeric-method
#' altExp<-,SingleCellExperiment,character-method
#' altExps<-,SingleCellExperiment-method
#' altExpNames<-,SingleCellExperiment,character-method
#' [,SummarizedExperimentByColumn,ANY,ANY,ANY-method
#' [<-,SummarizedExperimentByColumn,ANY,ANY,ANY-method
#' c,SummarizedExperimentByColumn-method
#' length,SummarizedExperimentByColumn-method
#' names,SummarizedExperimentByColumn-method
#' names<-,SummarizedExperimentByColumn-method
#' removeAltExps
#' mainExpName
#' mainExpName,SingleCellExperiment-method
#' mainExpName<-
#' mainExpName<-,SingleCellExperiment,character_OR_NULL-method
#'
#' % Dumping the SEBC methods here, so that check doesn't complain.
NULL

.alt_key <- "altExps"

#' @export
setMethod("altExpNames", "SingleCellExperiment", function(x) {
    .get_internal_names(x, 
        getfun=int_colData, 
        key=.alt_key)
})

#' @export
#' @importFrom S4Vectors endoapply
setMethod("altExps", "SingleCellExperiment", function(x, withDimnames=TRUE, withColData=FALSE) {
    value <- .get_internal_all(x, 
        getfun=int_colData, 
        key=.alt_key)

    value <- endoapply(value, .get_se)

    if (withDimnames || withColData) {
        for (i in seq_along(value)) {
            value[[i]] <- .fill_altexps_info(x, value[[i]], 
                withDimnames=withDimnames, withColData=withColData)
        }
    }

    value
})

#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment colData colData<-
.fill_altexps_info <- function(x, out, withDimnames, withColData) {
    if (withDimnames) {
        colnames(out) <- colnames(x)
    }

    if (withColData) {
        # We prepend to make it easier to reverse this operation in the setters.
        # Any new colData fields are added to the right side, so we just have
        # to check the harder-to-modify left side when reversing.
        prep <- cbind(colData(x), colData(out))
        rownames(prep) <- colnames(out)
        colData(out) <- prep
    }

    out
}

#' @export
setMethod("altExp", c("SingleCellExperiment", "missing"), function(x, e, withDimnames=TRUE, withColData=FALSE) {
    .get_internal_missing(x, 
        basefun=altExp, 
        namefun=altExpNames, 
        funstr="altExp",
        withDimnames=withDimnames,
        withColData=withColData)
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "numeric"), function(x, e, withDimnames=TRUE, withColData=FALSE) {
    out <- .get_internal_integer(x, e,
        getfun=int_colData,
        key=.alt_key,
        funstr="altExp",
        substr="e")

    out <- .get_se(out)

    .fill_altexps_info(x, out, withDimnames=withDimnames, withColData=withColData)
})

#' @export
setMethod("altExp", c("SingleCellExperiment", "character"), function(x, e, withDimnames=TRUE, withColData=FALSE) {
    out <- .get_internal_character(x, e,
        getfun=int_colData,
        key=.alt_key,
        funstr="altExp",
        substr="e",
        namestr="altExpNames")

    out <- .get_se(out)

    .fill_altexps_info(x, out, withDimnames=withDimnames, withColData=withColData)
})

#' @export
setReplaceMethod("altExpNames", c("SingleCellExperiment", "character"), function(x, value) {
    .set_internal_names(x, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.alt_key)
})

#' @importFrom SummarizedExperiment colData colData<-
.check_altexp_columns <- function(main, alt, withDimnames, withColData, fun='altExp', vname="value") {
    if (!is.null(alt)) {
        if (withDimnames) {
            if (!identical(colnames(main), colnames(alt))) {
                msg <- paste0("'colnames(", vname, ")' are not the same as 'colnames(x)' for '", 
                    fun, "<-'. This will be an error in the next release of Bioconductor.")
                warning(paste(strwrap(msg), collapse="\n"))
            }
        }

        if (withColData) {
            main.cd <- colData(main)
            ncd <- ncol(main.cd)
            alt.cd <- colData(alt)
            acd <- ncol(alt.cd)
            keep <- seq_len(acd) <= ncd

            if (acd < ncd || !identical(alt.cd[,keep,drop=FALSE], main.cd)) {
                warning("left-most columns of 'colData(", vname, ")' should be the same as 'colData(x)' when 'withColData=TRUE'")
            } else {
                colData(alt) <- alt.cd[,!keep,drop=FALSE]
            }
        }
    }

    alt
}

#' @export
#' @importClassesFrom S4Vectors SimpleList
setReplaceMethod("altExps", "SingleCellExperiment", function(x, withDimnames=TRUE, withColData=FALSE, ..., value) {
    if (withDimnames || withColData) {
        for (v in seq_along(value)) {
            value[[v]] <- .check_altexp_columns(x, value[[v]], 
                withDimnames=withDimnames, withColData=withColData, 
                fun='altExps', vname=sprintf("value[[%s]]", v))
        }
    }

    .set_internal_all(x, value, 
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.alt_key,
        convertfun=SummarizedExperimentByColumn,
        xdimfun=ncol,
        vdimfun=length,
        funstr="altExps",
        xdimstr="ncol",
        vdimstr="columns")
})

#' @export
removeAltExps <- function(x) {
    altExps(x) <- NULL
    x
}

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "missing"), function(x, e, withDimnames=TRUE, withColData=FALSE, ..., value) {
    .set_internal_missing(x, value,
        withDimnames=withDimnames,
        withColData=withColData,
        basefun=`altExp<-`,
        namefun=altExpNames)
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "numeric"), function(x, e, withDimnames=TRUE, withColData=FALSE, ..., value) {
    value <- .check_altexp_columns(x, value, withDimnames=withDimnames, withColData=withColData)

    .set_internal_numeric(x, e, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.alt_key,
        convertfun=SummarizedExperimentByColumn,
        xdimfun=ncol,
        vdimfun=length,
        funstr="altExp",
        xdimstr="ncol",
        vdimstr="columns",
        substr="e")
})

#' @export
setReplaceMethod("altExp", c("SingleCellExperiment", "character"), function(x, e, withDimnames=TRUE, withColData=FALSE, ..., value) {
    value <- .check_altexp_columns(x, value, withDimnames=withDimnames, withColData=withColData)

    .set_internal_character(x, e, value, 
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.alt_key,
        convertfun=SummarizedExperimentByColumn,
        xdimfun=ncol, 
        vdimfun=length,
        funstr="altExp", 
        xdimstr="ncol",
        vdimstr="columns", 
        substr="e") 
})

#' @export
setMethod("mainExpName", "SingleCellExperiment", function(x) {
    int_metadata(x)$mainExpName
})

#' @export
setReplaceMethod("mainExpName", c("SingleCellExperiment", "character_OR_NULL"), function(x, value) {
    int_metadata(x)$mainExpName <- value
    x
})
