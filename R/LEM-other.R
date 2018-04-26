#############################################
# This file defines the non-get-set methods for the LinearEmbeddingMatrix class.
#############################################

#############################################
# Sets the validity checker.

.le_validity <- function(object) {
    msg <- NULL

    # Check dimensions
    sf <- sampleFactors(object, withDimnames=FALSE)
    fl <- featureLoadings(object, withDimnames=FALSE)
    fd <- factorData(object, withDimnames=FALSE)

    if(NCOL(sf) != NCOL(fl)) {
        msg <- c(msg, "'sampleFactors' and 'featureLoadings' must have the same number of columns")
    }

    if(NROW(fd) != NCOL(sf)) {
        msg <- c(msg, "'factorData' must have one row per factor")
    }

    rn <- rownames(object)
    if (!is.null(rn) && length(rn) != NROW(sf)) {
        msg <- c(msg, "'length(NAMES)' must be NULL or equal to 'nrow(object)'")
    }
    
    if (length(msg)) { return(msg) }
    return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("LinearEmbeddingMatrix", .le_validity)

#############################################
# Sets the show method.

.le_show <- function(object) {

    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    ## metadata
    expt <- names(metadata(object))
    if (is.null(expt)) {
        expt <- character(length(metadata(object)))
    }
    scat("metadata(%d): %s\n", expt)

    ## rownames
    rownames <- rownames(object)
    if(is.null(rownames)) {
        cat("rownames: NULL\n")
    } else {
        scat("rownames(%d): %s\n", rownames(object))
    }

    ## colnames
    colnames <- colnames(object)
    if(is.null(colnames)) {
        cat("colnames: NULL\n")
    } else {
        scat("colnames(%d): %s\n", colnames(object))
    }

    ## factorData
    scat("factorData names(%d): %s\n", names(factorData(object)))
}

#' @export
setMethod("show", "LinearEmbeddingMatrix", .le_show)

#############################################
# Defines a constructor

#' @export
#' @importClassesFrom S4Vectors DataFrame
LinearEmbeddingMatrix <- function(sampleFactors = matrix(nrow = 0, ncol = 0),
                            featureLoadings = matrix(nrow = 0, ncol = 0),
                            factorData = NULL,
                            metadata = list()) {
    if (is.null(factorData)) {
        factorData <- new("DataFrame", nrows = ncol(sampleFactors))
    }
    out <- new("LinearEmbeddingMatrix",
               sampleFactors = sampleFactors,
               featureLoadings = featureLoadings,
               factorData = factorData,
               metadata = as.list(metadata)
               )
    return(out)
}

#############################################
# Define subsetting methods.

#' @export
setMethod("[", c("LinearEmbeddingMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    temp_sf <- sampleFactors(x, withDimnames=FALSE)
    temp_fl <- featureLoadings(x, withDimnames=FALSE)
    temp_fd <- factorData(x, withDimnames=FALSE)
    temp_rn <- rownames(x)

    if(!missing(i)) {
        i <- .convert_subset_index(i, rownames(x))
        temp_sf <- temp_sf[i,,drop=drop]
        temp_rn <- temp_rn[i]
    }

    if(!missing(j)) {
        j <- .convert_subset_index(j, colnames(x))
        temp_sf <- temp_sf[,j,drop=drop]
        temp_fl <- temp_fl[,j,drop=drop]
        temp_fd <- temp_fd[j,,drop=FALSE]
    }

    # Returning a vector, a la drop=TRUE for a full matrix.
    if (is.null(dim(temp_sf))) {
        return(temp_sf)
    }

    BiocGenerics:::replaceSlots(x, sampleFactors = temp_sf,
               featureLoadings = temp_fl, 
               factorData = temp_fd, 
               NAMES = temp_rn, check=FALSE)
})

#' @export
setMethod("[<-", c("LinearEmbeddingMatrix", "ANY", "ANY", "LinearEmbeddingMatrix"), function(x, i, j, ..., value) {
    temp_sf <- sampleFactors(x, withDimnames=FALSE)
    temp_fl <- featureLoadings(x, withDimnames=FALSE)
    temp_fd <- factorData(x, withDimnames=FALSE)

    # i is samples
    # j is factors
    if (!missing(i)) {
        i <- .convert_subset_index(i, rownames(x))
    }
    if (!missing(j)) {
        j <- .convert_subset_index(j, colnames(x))
    }

    # Inserting sample factors.
    if (missing(i) && missing(j)) {
        temp_sf <- sampleFactors(value, withDimnames=FALSE)
    } else if (missing(i)) {
        temp_sf[,j] <- sampleFactors(value, withDimnames=FALSE)
    } else if (missing(j)) {
        temp_sf[i,] <- sampleFactors(value, withDimnames=FALSE)
    } else {
        temp_sf[i,j] <- sampleFactors(value, withDimnames=FALSE)
    }

    # Dealing with the factorData, featureLoadings.
    if (missing(i) && missing(j)) {
        temp_fl <- featureLoadings(value)
        temp_fd <- factorData(value)
    } else if (!missing(j)) {
        temp_fl[,j] <- featureLoadings(value)
        temp_fd[j,] <- factorData(value)
    }
   
    # Dealing with NAMES. 
    temp_rn <- rownames(x)
    if (!missing(i)) {
        temp_rn[i] <- rownames(value)
    }

    BiocGenerics:::replaceSlots(x, sampleFactors = temp_sf,
               featureLoadings = temp_fl, 
               factorData = temp_fd, 
               NAMES = temp_rn, check=FALSE)
})

#############################################
# Defining the combining methods.

#' @importFrom BiocGenerics rbind
setMethod("rbind", "LinearEmbeddingMatrix", function(..., deparse.level=1) {
    args <- list(...)
    x <- args[[1]]
    all_sf <- lapply(args, sampleFactors)
    all_sf <- do.call(rbind, all_sf)

    # Checking what to do with names.
    all_rn <- lapply(args, rownames)
    unnamed <- sapply(all.names, is.null)
    for (u in which(unnamed)) {
        all_rn[[u]] <- character(length(args[[u]]))
    }
    all_rn <- unlist(all_rn)

    # Replacing the relevant slots.
    BiocGenerics:::replaceSlots(x, NAMES=all_rn, sampleFactors=all_sf, check=FALSE)
})

#' @importFrom BiocGenerics cbind
setMethod("cbind", "LinearEmbeddingMatrix", function(..., deparse.level=1) {
    args <- list(...)
    x <- args[[1]]

    all_sf <- lapply(args, sampleFactors, withDimnames=FALSE)
    all_fd <- lapply(args, factorData, withDimnames=FALSE)
    all_fl <- lapply(args, featureLoadings, withDimnames=FALSE)

    BiocGenerics:::replaceSlots(x, sampleFactors = do.call(cbind, all_sf),
               featureLoadings = do.call(cbind, all_fl),
               factorData = do.call(rbind, all_fd), check=FALSE)
})

