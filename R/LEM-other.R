#############################################
# This file defines the non-get-set methods for the LinearEmbeddingMatrix class.
#############################################

#############################################
# Sets the validity checker.

.le_validity <- function(object) {
    msg <- NULL

    # Check dimensions
    sf <- sampleFactors(object)
    fl <- featureLoadings(object)
    fd <- factorData(object)

    if(NCOL(sf) != NCOL(fl)) {
        msg <- c(msg, "'sampleFactors' and 'featureLoadings' must have the same number of columns")
    }

    if(NROW(fd) != NCOL(sf)) {
        msg <- c(msg, "'factorData' must have one row per factor")
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("LinearEmbeddingMatrix", .le_validity)

#############################################
# Sets the show method.

.le_show <- function(object) {
    cat("class: LinearEmbeddingMatrix", "\n")
    cat(sprintf("Number of factors: %d\n", ncol(sampleFactors(object))))
    cat(sprintf("Number of samples: %d\n", nrow(sampleFactors(object))))
    cat(sprintf("Number of features: %d\n", nrow(featureLoadings(object))))
}

#' @export
setMethod("show", "LinearEmbeddingMatrix", .le_show)

#############################################
# Defines a constructor

#' @export
#' @importClassesFrom S4Vectors DataFrame
LinearEmbeddingMatrix <- function(sampleFactors = matrix(nrow = 0, ncol = 0),
                            featureLoadings = matrix(nrow = 0, ncol = 0),
                            factorData = NULL) {
    if (is.null(factorData)) {
        factorData <- new("DataFrame", nrows = ncol(sampleFactors))
    }
    out <- new("LinearEmbeddingMatrix",
               sampleFactors = sampleFactors,
               featureLoadings = featureLoadings,
               factorData = factorData
               )
    return(out)
}

#############################################
# Define subsetting methods.

#' @export
setMethod("[", c("LinearEmbeddingMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    temp_sf <- sampleFactors(x)
    temp_fl <- featureLoadings(x)
    temp_fd <- factorData(x)

    if(!missing(i)) {
        temp_sf <- temp_sf[i,,drop=drop]
    }

    if(!missing(j)) {
        temp_sf <- temp_sf[,j,drop=drop]
        temp_fl <- temp_fl[,j,drop=drop]
        temp_fd <- temp_fd[j,,drop=FALSE]
    }

    # Returning a vector, a la drop=TRUE for a full matrix.
    if (is.null(dim(temp_sf))) {
        return(temp_sf)
    }

    initialize(x, sampleFactors = temp_sf,
               featureLoadings = temp_fl,
               factorData = temp_fd)
})

#' @export
setMethod("[<-", c("LinearEmbeddingMatrix", "ANY", "ANY", "LinearEmbeddingMatrix"), function(x, i, j, ..., value) {
    temp_sf <- sampleFactors(x)
    temp_fl <- featureLoadings(x)
    temp_fd <- factorData(x)

    # i is samples
    # j is factors
    if (missing(i) && missing(j)) {
        temp_sf <- sampleFactors(value)
        temp_fl <- featureLoadings(value)
        temp_fd <- factorData(value)
    } else if (missing(i)) {
        temp_sf[,j] <- sampleFactors(value)
    } else if (missing(j)) {
        temp_sf[i,] <- sampleFactors(value)
    } else {
        temp_sf[i,j] <- sampleFactors(value)
    }

    if (!missing(j)) {
        temp_fl[,j] <- featureLoadings(value)
        temp_fd[j,] <- factorData(value)
    }

    initialize(x, sampleFactors = temp_sf,
               featureLoadings = temp_fl,
               factorData = temp_fd)
})

#############################################
# Defining the combining methods.

#' @importFrom BiocGenerics rbind
setMethod("rbind", "LinearEmbeddingMatrix", function(..., deparse.level=1) {
    args <- list(...)
    ans <- args[[1]]
    all_sf <- lapply(args, sampleFactors)
    sampleFactors(ans) <- do.call(rbind, all_sf)
    return(ans)
})

#' @importFrom BiocGenerics cbind
setMethod("cbind", "LinearEmbeddingMatrix", function(..., deparse.level=1) {
    args <- list(...)
    ans <- args[[1]]
    args <- args[-1]

    all_sf <- all_fd <- all_fl <- vector("list", length(args)+1L)
    all_sf[[1]] <- sampleFactors(ans)
    all_fd[[1]] <- factorData(ans)
    all_fl[[1]] <- featureLoadings(ans)

    for (x in seq_along(args)) {
        current <- args[[x]]
        all_sf[[x+1]] <- sampleFactors(current)
        all_fd[[x+1]] <- factorData(current)
        all_fl[[x+1]] <- featureLoadings(current)
    }

    initialize(ans, sampleFactors = do.call(cbind, all_sf),
               featureLoadings = do.call(cbind, all_fl),
               factorData = do.call(rbind, all_fd))
})

