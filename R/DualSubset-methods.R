# Basic methods for the DualSubset class,
# an internal class that powers the *Pairs methods.

.get_hits <- function(x) x@hits

DualSubset <- function(hits) new("DualSubset", hits=sort(hits))

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
    SelfHits(mat@i + 1L, mat@j + 1L, nnode=nrow(mat), x=mat@x)
}

#' @importFrom S4Vectors mcols mcols<- findMatches queryHits subjectHits SelfHits
setMethod("[", "DualSubset", function(x, i, j, ..., drop=FALSE) {
    p <- .get_hits(x)
    i <- normalizeSingleBracketSubscript(i, x)

    mq <- findMatches(i, queryHits(p))
    ms <- findMatches(i, subjectHits(p))
    mcom <- findMatches(subjectHits(mq), subjectHits(ms))

    left <- queryHits(mq)[queryHits(mcom)]
    right <- queryHits(ms)[subjectHits(mcom)]
    index <- subjectHits(mq)[queryHits(mcom)] # same as subjectHits(ms)[subjectHits(mcom)]

    o <- order(left, right)
    hits2 <- SelfHits(left[o], right[o], nnode=length(i))
    mcols(hits2) <- mcols(p)[index[o],,drop=FALSE]
    initialize(x, hits=hits2)
})

#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits SelfHits
setReplaceMethod("[", "DualSubset", function(x, i, j, ..., value) {
    p <- .get_hits(x)
    i <- normalizeSingleBracketSubscript(i, x)

    # Filtering out all elements to be replaced.
    p <- p[!(queryHits(p) %in% i & subjectHits(p) %in% i)]

    pv <- .get_hits(value)
    new.q <- i[queryHits(pv)]
    new.s <- i[subjectHits(pv)]

    total.p <- c(queryHits(p), new.q)
    total.s <- c(subjectHits(p), new.s)
    total.m <- rbind(mcols(p), mcols(pv))

    o <- order(total.p, total.s)
    hits2 <- SelfHits(total.p[o], total.s[o], nnode=length(x))
    mcols(hits2) <- total.m[o,,drop=FALSE]
    initialize(x, hits=hits2)
})

#' @importFrom utils tail 
#' @importFrom S4Vectors queryHits subjectHits SelfHits mcols mcols<- 
setMethod("c", "DualSubset", function(x, ...) {
    everything <- list(x, ...)
    shift <- 0L

    all.first <- all.second <- all.values <- vector("list", length(everything))
    for (i in seq_along(everything)) {
        current <- .get_hits(everything[[i]])
        contribution <- nnode(current)
        all.first[[i]] <- queryHits(current) + shift
        all.second[[i]] <- subjectHits(current) + shift
        all.values[[i]] <- mcols(current)
        shift <- shift + contribution
    }

    final <- SelfHits(unlist(all.first), unlist(all.second), nnode=shift)
    mcols(final) <- do.call(rbind, all.values)
    initialize(x, hits=final)
})
