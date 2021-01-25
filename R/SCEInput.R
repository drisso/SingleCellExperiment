#' The SCEInput class
#'
#' The SCEInput classes provide a flexible approach to extracting different types of data from a \linkS4class{SingleCellExperiment} object.
#' Each subclass specifies a different piece of data to fetch from the SingleCellExperiment in a generic manner,
#' as well as storing custom arguments that can be to process that particular piece of data.
#' This is helpful for general-purpose functions like \code{\link{applySCE}},
#' as well as complex functions that use multiple pieces of information from a SingleCellExperiment object.
#'
#' @section Constructors:
#' One constructor is available for each SCEInput subclass:
#' \itemize{
#' \item \code{MainExpInput(...)}
#' \item \code{AssayInput(assay, ...)}
#' \item \code{ReducedDimInput(type, ...)}
#' \item \code{AltExpInput(experiment, ...)}
#' \item \code{AltAssayInput(experiment, assay, ...)}
#' \item \code{AltReducedDimInput(experiment, type, ...)}
#' }
#' 
#' With the following arguments:
#' \describe{
#' \item{...}{Aguments (usually named) to pass to other functions, e.g., \code{FUN} in \code{\link{applySCE}}.}
#' \item{assay}{String or integer specifying the assay to use.}
#' \item{type}{String or integer specifying the \code{\link{reducedDim}} to use.}
#' \item{experiment}{String or integer specifying the \code{\link{altExp}} to use.}
#' }
#'
#' @section Getting the data:
#' \code{getInput(x, input)} will return the specified piece of data from a \linkS4class{SingleCellExperiment} \code{x},
#' based on the SCEInput object \code{input}.
#' \itemize{
#' \item If \code{input} is a MainExpInput, \code{x} is returned directly.
#' \item If \code{input} is a AssayInput, the specified assay of \code{x} is returned. 
#' \item If \code{input} is a ReducedDimInput, the specified entry of \code{\link{reducedDims}(x)} is returned. 
#' \item If \code{input} is an AltExpInput, the specified alternative Experiment in \code{x} is returned.
#' \item If \code{input} is an AltAssayInput, the specified assay of the specified alternativeExperiment is returned.
#' \item If \code{input} is an AltReducedDimInput, the specified \code{\link{reducedDims}} of the specified alternativeExperiment is returned.
#' }
#'
#' @section Input properties:
#' In the following code snippets, \code{input} is a SCEInput object.
#'
#' \code{inputArguments(input)} will return a list of arguments from \code{input}.
#'
#' \code{isReducedDimInput(input)} will return \code{TRUE} for ReducedDimInput and AltExpReducedInput objects, and \code{FALSE} otherwise.
#'
#' \code{isInputSameShape(input)} will return \code{TRUE} for all \code{input} where \code{getInput(x, input)} is guaranteed to have the same dimensions as \code{x}.
#'
#' \code{isInputAltExp(input)} will return \code{TRUE} for all \code{input} where \code{getInput(x, input)} pulls data from \code{altExps(x)}.
#'
#' @section Constructing lists of SCEInputs:
#' In the following code snippets, \code{x} is a SingleCellExperiment object.
#' 
#' \code{makeAllExpInputs(x, include.main=TRUE)} will return a list containing AltExpInputs for all \code{\link{altExpNames}(x)}.
#' If \code{include.main=TRUE}, a MainExpInput is prepended to the list.
#'
#' \code{makeSameAssayInputs(x, assay, include.main=TRUE)} will return a list containing AltAssayInputs with the specified \code{assay} for all \code{\link{altExpNames}(x)}.
#' If \code{include.main=TRUE}, an AssayInput is prepended to the list.
#' 
#' 
#' @author Aaron Lun
#'
#' @examples
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' pca <- matrix(runif(ncells*5), ncells)
#' sce <- SingleCellExperiment(assays=list(counts=u),
#'     reducedDims=SimpleList(PCA=pca))
#' sce.copy <- SingleCellExperiment(assays=list(counts=u*10),
#'     reducedDims=SimpleList(PCA=pca/10))
#' altExp(sce, "BLAH") <- sce.copy
#'
#' getInput(sce, MainExpInput())
#' str(getInput(sce, AssayInput("counts")))
#' str(getInput(sce, ReducedDimInput("PCA")))
#' getInput(sce, AltExpInput("BLAH"))
#' str(getInput(sce, AltAssayInput("BLAH", "counts")))
#' str(getInput(sce, AltReducedDimInput("BLAH", "PCA")))
#' 
#' @seealso
#' \code{\link{applySCE}}, for an example of how this can be used in a generic manner.
#'
#' @docType class
#' @name SCEInput
#' @aliases
#' SCEInput
#' MainExpInput
#' AssayInput
#' ReducedDimInput
#' AltExpInput
#' AltAssayInput
#' AltReducedDimInput
#' SCEInput-class
#' MainExpInput-class
#' AssayInput-class
#' ReducedDimInput-class
#' AltExpInput-class
#' AltAssayInput-class
#' AltReducedDimInput-class
#' show
#' show,SCEInput-method
#' getInput
#' getInput,SingleCellExperiment,MainExpInput-method
#' getInput,SingleCellExperiment,AssayInput-method
#' getInput,SingleCellExperiment,ReducedDimInput-method
#' getInput,SingleCellExperiment,AltExpInput-method
#' getInput,SingleCellExperiment,AltAssayInput-method
#' getInput,SingleCellExperiment,AltReducedDimInput-method
#' isInputReducedDim
#' isInputReducedDim,SCEInput-method
#' isInputReducedDim,ReducedDimInput-method
#' isInputReducedDim,AltReducedDimInput-method
#' isInputSameShape
#' isInputSameShape,SCEInput-method
#' isInputSameShape,ReducedDimInput-method
#' isInputSameShape,AltReducedDimInput-method
#' isInputAltExp
#' isInputAltExp,SCEInput-method
#' isInputAltExp,AltExpInput-method
#' inputArguments
#' makeAllExpInputs 
#' makeSameAssayInputs
NULL

#' @export
MainExpInput <- function(...) {
    new("MainExpInput", arguments=list(...))
}

#' @export
AssayInput <- function(assay, ...) {
    new("AssayInput", assay=assay, arguments=list(...))
}

.check_valid_index <- function(value, name) {
    if (is.character(value) || is.numeric(value)) {
        if (length(value)==1 && !is.na(value)) {
            return(TRUE)
        }
        sprintf("'%s' must be a non-NA scalar", name)
    } else {
        sprintf("'%s' must be a string or integer", name)
    }
}

setValidity2("AssayInput", function(object){
    .check_valid_index(object@assay, "assay")
})

#' @export
ReducedDimInput <- function(type, ...) {
    new("ReducedDimInput", type=type, arguments=list(...))
}

setValidity2("ReducedDimInput", function(object){
    .check_valid_index(object@type, "type")
})

#' @export
AltExpInput <- function(experiment, ...) {
    new("AltExpInput", experiment=experiment, arguments=list(...))
}

setValidity2("AltExpInput", function(object){
    .check_valid_index(object@experiment, "experiment")
})

#' @export
AltAssayInput <- function(experiment, assay, ...) {
    new("AltAssayInput", experiment=experiment, assay=assay, arguments=list(...))
}

setValidity2("AltAssayInput", function(object){
    .check_valid_index(object@assay, "assay")
})

#' @export
AltReducedDimInput <- function(experiment, type, ...) {
    new("AltReducedDimInput", experiment=experiment, type=type, arguments=list(...))
}

setValidity2("AltReducedDimInput", function(object){
    .check_valid_index(object@type, "type")
})

###################################################

#' @export
setMethod("show", "SCEInput", function(object) {
    cat(paste0("class: ", class(object), "\n"))
    .input_subshow(object)
    coolcat("arguments(%d): %s\n", names(inputArguments(object)))
})

setGeneric(".input_subshow", function(object) standardGeneric(".input_subshow"))

setMethod(".input_subshow", "MainExpInput", function(object) NULL)

setMethod(".input_subshow", "AssayInput", function(object) {
    cat(paste0("assay: ", object@assay, "\n")) 
})

setMethod(".input_subshow", "ReducedDimInput", function(object) {
    cat(paste0("type: ", object@type, "\n")) 
})

setMethod(".input_subshow", "AltExpInput", function(object) {
    cat(paste0("experiment: ", object@experiment, "\n")) 
})

setMethod(".input_subshow", "AltAssayInput", function(object) {
    callNextMethod()
    cat(paste0("assay: ", object@assay, "\n")) 
})

setMethod(".input_subshow", "AltReducedDimInput", function(object) {
    callNextMethod()
    cat(paste0("type: ", object@type, "\n")) 
})

###################################################

#' @export
setMethod("getInput", c(x="SingleCellExperiment", input="MainExpInput"), function(x, input) {
    x 
})

#' @export
setMethod("getInput", c(x="SingleCellExperiment", input="AssayInput"), function(x, input) {
    assay(x, input@assay)  
})

#' @export
setMethod("getInput", c(x="SingleCellExperiment", input="ReducedDimInput"), function(x, input) {
    reducedDim(x, input@type)  
})

#' @export
setMethod("getInput", c(x="SingleCellExperiment", input="AltExpInput"), function(x, input) {
    altExp(x, input@experiment)
})

#' @export
setMethod("getInput", c(x="SingleCellExperiment", input="AltAssayInput"), function(x, input) {
    assay(altExp(x, input@experiment), input@assay)  
})

#' @export
setMethod("getInput", c(x="SingleCellExperiment", input="AltReducedDimInput"), function(x, input) {
    reducedDim(altExp(x, input@experiment), input@type)  
})

###################################################

#' @export
setMethod("isInputReducedDim", "SCEInput", function(input) {
    FALSE 
})

#' @export
setMethod("isInputReducedDim", "ReducedDimInput", function(input) {
    TRUE
})

#' @export
setMethod("isInputReducedDim", "AltReducedDimInput", function(input) {
    TRUE
})

#' @export
setMethod("isInputSameShape", "SCEInput", function(input) {
    TRUE 
})

#' @export
setMethod("isInputSameShape", "ReducedDimInput", function(input) {
    FALSE 
})

#' @export
setMethod("isInputSameShape", "AltReducedDimInput", function(input) {
    FALSE
})

#' @export
setMethod("isInputAltExp", "SCEInput", function(input) {
    FALSE
})

#' @export
setMethod("isInputAltExp", "AltExpInput", function(input) {
    TRUE
})

#' @export
inputArguments <- function(input) {
    input@arguments
}

###################################################

#' @export
makeAllExpInputs <- function(x, include.main=TRUE) {
    output <- lapply(altExpNames(x), AltExpInput)
    if (include.main) {
        output <- c(list(MainExpInput()), output)
    }
    output
}

#' @export
makeSameAssayInputs <- function(x, assay, include.main=TRUE) {
    output <- lapply(altExpNames(x), AltAssayInput, assay=assay)
    if (include.main) {
        output <- c(list(AssayInput(assay=assay)), output)
    }
    output
}
