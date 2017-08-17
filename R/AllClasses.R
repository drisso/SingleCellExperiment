# Defines the SingleCellExperiment class.

setClass("SingleCellExperiment",
         slots=c(int_elementMetadata = "DataFrame",
                 int_colData = "DataFrame",
                 int_metadata = "list",
                 reducedDims = "SimpleList"),
         contains = "RangedSummarizedExperiment")

#############################################
# Sets the validity checker.

.sce_validity <- function(object) {
    msg <- NULL

    # Checking dimensions of reduced coordinates.
    rd <- reducedDims(object)
    if (length(rd)) {
        if (any(unlist(lapply(rd, .not_reddim_mat, object=object)))) {
            msg <- c(msg,
                     "each element of 'reducedDims' must be a matrix with nrow equal to 'ncol(object)'")
        }
    }

    # Checking dimensions of internal objects.
    if (nrow(int_elementMetadata(object))!=nrow(object)) {
        msg <- c(msg, "'nrow' of internal 'rowData' not equal to 'nrow(object)'")
    }
    if (nrow(int_colData(object))!=ncol(object)) {
        msg <- c(msg, "'nrow' of internal 'colData' not equal to 'ncol(object)'")
    }

    # Checking spike-in names are present and accounted for.
    spike.fields <- .get_spike_field(spikeNames(object), check=FALSE)
    lost.spikes <- ! spike.fields %in% colnames(int_elementMetadata(object))
    if (any(lost.spikes)) {
        was.lost <- spikeNames(object)[lost.spikes][1]
        msg <- c(msg, sprintf("no field specifying rows belonging to spike-in set '%s'", was.lost))
    }

    # Checking version.
    if (objectVersion(object) < "0.98.0") {
        msg <- c(msg, "object is out of date, update with 'updateSCE(object)'")
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

.not_reddim_mat <- function(val, object) {
    return(!is.matrix(val) || nrow(val)!=ncol(object));
}

setValidity2("SingleCellExperiment", .sce_validity)

#############################################
# Sets the show method.

scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

.sce_show <- function(object) {
    callNextMethod()
    scat("reduced(%d): %s\n", names(reducedDims(object)))
    scat("spikes(%d): %s\n", spikeNames(object))
}

setMethod("show", "SingleCellExperiment", .sce_show)

#############################################
# Defines a constructor.

SingleCellExperiment <- function(..., reducedDims=SimpleList()) {
    se <- SummarizedExperiment(...)
    if(!is(se, "RangedSummarizedExperiment")) {
      rse <- as(se, "RangedSummarizedExperiment")
      rowData(rse) <- rowData(se)
    } else {
      rse <- se
    }
    out <- new("SingleCellExperiment", rse, reducedDims=SimpleList(),
               int_elementMetadata=DataFrame(matrix(0, nrow(se), 0)),
               int_colData=DataFrame(matrix(0, ncol(se), 0)),
               int_metadata=list(version=packageVersion("SingleCellExperiment")))
    reducedDims(out) <- reducedDims
    return(out)
}

#############################################
# Define subsetting methods.

setMethod("[", c("SingleCellExperiment", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) {
        ii <- .convert_subset_index(i, rownames(x))
        int_elementMetadata(x) <- int_elementMetadata(x)[ii,,drop=FALSE]
    }

    if (!missing(j)) {
        jj <- .convert_subset_index(j, colnames(x))
        int_colData(x) <- int_colData(x)[jj,,drop=FALSE]
        rd <- reducedDims(x)
        for (mode in seq_along(rd)) { rd[[mode]] <- rd[[mode]][jj,,drop=FALSE] }
        int_reducedDims(x) <- rd
    }

    callNextMethod()
})

setMethod("[<-", c("SingleCellExperiment", "ANY", "ANY", "SingleCellExperiment"), function(x, i, j, ..., value) {
    if (missing(i) && missing(j)) {
        int_elementMetadata(x) <- int_elementMetadata(value)
        int_colData(x) <- int_colData(value)
        int_reducedDims(x) <- reducedDims(value)
    }

    if (!missing(i)) {
        ii <- .convert_subset_index(i, rownames(x))
        sout <- .standardize_DataFrames(first=int_elementMetadata(x), last=int_elementMetadata(value))
        sout$first[ii,] <- sout$last
        int_elementMetadata(x) <- sout$first
    }

    if (!missing(j)) {
        jj <- .convert_subset_index(j, colnames(x))
        sout <- .standardize_DataFrames(first=int_colData(x), last=int_colData(value))
        sout$first[jj,] <- sout$last
        int_colData(x) <- sout$first

        rdout <- .standardize_reducedDims(first=x, last=value)
        rd <- rdout$first
        rdv <- rdout$last
        for (mode in seq_along(rd)) {
            rd[[mode]][jj,] <- rdv[[mode]]
        }
        int_reducedDims(x) <- rd
    }

    int_metadata(x) <- int_metadata(value)
    callNextMethod()
})

setMethod("subset", "SingleCellExperiment", function(x, i, j) {
    x[i, j]
})

#############################################
# Defining the combining methods.

setMethod("cbind", "SingleCellExperiment", function(..., deparse.level=1) {
    args <- unname(list(...))
    base <- do.call(cbind, lapply(args, function(x) { as(x, "RangedSummarizedExperiment") }))

    all.col.data <- lapply(args, int_colData)
    sout <- do.call(.standardize_DataFrames, all.col.data)
    new.col.data <- do.call(rbind, sout)

    all.rd <- do.call(.standardize_reducedDims, args)
    new.rd <- SimpleList(do.call(mapply, c(all.rd, FUN=rbind, SIMPLIFY=FALSE)))

    ans <- args[[1]]
    new("SingleCellExperiment", base, int_colData=new.col.data, int_elementMetadata=int_elementMetadata(ans),
        int_metadata=int_metadata(ans), reducedDims=new.rd)
})

setMethod("rbind", "SingleCellExperiment", function(..., deparse.level=1) {
    args <- unname(list(...))
    base <- do.call(rbind, lapply(args, function(x) { as(x, "RangedSummarizedExperiment") }))

    all.row.data <- lapply(args, int_elementMetadata)
    sout <- do.call(.standardize_DataFrames, all.row.data)
    new.row.data <- do.call(rbind, sout)

    ans <- args[[1]]
    new("SingleCellExperiment", base, int_colData=int_colData(ans), int_elementMetadata=new.row.data,
        int_metadata=int_metadata(ans), reducedDims=reducedDims(ans))
})

