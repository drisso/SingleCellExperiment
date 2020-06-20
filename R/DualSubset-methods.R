# Basic methods for the DualSubset class,
# an internal class that powers the *Pairs methods.

.get_hits <- function(x) x@hits

DualSubset <- function(hits, values) new("DualSubset", hits=hits, values=values)

#' @importFrom S4Vectors nnode
setMethod("length", "DualSubset", function(x) nnode(.get_hits(x)))

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors first second 
.hits2mat <- function(p) {
    sparseMatrix(i=first(p), j=second(p), x=seq_along(p), giveCsparse=FALSE)
}

#' @importFrom S4Vectors mcols SelfHits mcols<-
.mat2hits <- function(mat) {
    out <- SelfHits(mat@i, mat@j, nnode=nrow(mat))
    mcols(out)$value <- mat@x
    out
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
#' @importFrom S4Vectors first second SelfHits mcols mcols<-
setMethod("c", "DualSubset", function(x, ...) {
    everything <- list(x, ...)
    N <- sum(lengths(everything))
    shift <- 0L

    all.first <- all.second <- all.values <- vector("list", length(gathered))
    for (i in seq_along(gathered)) {
        current <- .get_hits(gathered[[i]])
        contribution <- nnode(current)
        all.first[[i]] <- first(current) + shift
        all.second[[i]] <- second(current) + shift
        all.values[[i]] <- mcols(current)$value
        shift <- shift + contribution
    }

    final <- SelfHits(unlist(all.first), unlist(all.second), nnode=N)
    mcols(final)$value <- unlist(all.values)
    initialize(x, hits=final)
})