# Basic methods for the DualSubset class,
# an internal class that powers the *Pairs methods.

.get_hits <- function(x) x@hits

DualSubset <- function(hits) new("DualSubset", hits=hits)

#' @importFrom S4Vectors nnode
setMethod("length", "DualSubset", function(x) nnode(.get_hits(x)))

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors queryHits subjectHits
.hits2mat <- function(p, x=seq_along(p)) {
    if (!is.logical(x) && !is.numeric(x) && !is.complex(x)) {
        stop("values of type '", typeof(x), "' are not supported in sparse matrices")
    }
    sparseMatrix(i=queryHits(p), j=subjectHits(p), x=x, 
        dims=rep(nnode(p), 2L), giveCsparse=FALSE)
}

#' @importClassesFrom Matrix dgTMatrix
#' @importFrom S4Vectors SelfHits 
.mat2hits <- function(mat) {
    mat <- as(mat, "dgTMatrix")
    SelfHits(mat@i + 1L, mat@j + 1L, nnode=nrow(mat), value=mat@x)
}

#' @importFrom S4Vectors mcols mcols<-
setMethod("[", "DualSubset", function(x, i, j, ..., drop=FALSE) {
    p <- .get_hits(x)
    mat <- .hits2mat(p)
    mat <- mat[i,i,drop=FALSE]
    p2 <- .mat2hits(mat)
    mcols(p2)$value <- mcols(p)$value[mcols(p2)$value]
    initialize(x, hits=p2)
})

#' @importFrom S4Vectors mcols mcols<-
setReplaceMethod("[", "DualSubset", function(x, i, j, ..., value) {
    p <- .get_hits(x)
    mat <- .hits2mat(p)
    pv <- .get_hits(value)
    matv <- -.hits2mat(pv)

    # Wiping out all pairs involving the elements to be replaced.
    mat[i,] <- 0
    mat[,i] <- 0
    mat[i,i] <- matv

    p2 <- .mat2hits(mat)
    index <- mcols(p2)$value
    use.left <- index > 0

    if (any(use.left)) {
        store <- mcols(p)$value[ifelse(use.left, index, 1)]
        store[!use.left] <- mcols(pv)$value[-index[!use.left]]
    } else if (any(!use.left)) {
        store <- mcols(pv)$value[-index]
    }

    mcols(p2)$value <- store
    initialize(x, hits=p2)
})

#' @importFrom utils tail 
#' @importFrom S4Vectors queryHits subjectHits SelfHits 
setMethod("c", "DualSubset", function(x, ...) {
    everything <- list(x, ...)
    shift <- 0L

    all.first <- all.second <- all.values <- vector("list", length(everything))
    for (i in seq_along(everything)) {
        current <- .get_hits(everything[[i]])
        contribution <- nnode(current)
        all.first[[i]] <- queryHits(current) + shift
        all.second[[i]] <- subjectHits(current) + shift
        all.values[[i]] <- mcols(current)$value
        shift <- shift + contribution
    }

    final <- SelfHits(unlist(all.first), unlist(all.second), nnode=shift, 
        value=unlist(all.values))
    initialize(x, hits=final)
})
