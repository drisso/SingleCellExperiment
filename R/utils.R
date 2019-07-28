.convert_subset_index <- function(subset, names) {
    if (is.character(subset)) {
        fmt <- "index out of bounds: %s"
        subset <- SummarizedExperiment:::.SummarizedExperiment.charbound(subset, names, fmt)
    }
    return(as.vector(subset))
}

.fill_int_fields <- function(df, fields) {
    for (field in fields) {
        if (!field %in% colnames(df)) {
            df[[field]] <- df[,0]
        }
    }
    df
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.create_shell_coldata <- function(x) {
    df <- int_colData(x)
    df <- .fill_int_fields(df, c(.red_key, .alt_key))
    SummarizedExperiment(colData=df)
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.create_shell_rowdata <- function(x) {
    SummarizedExperiment(rowData=int_elementMetadata(x))
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.create_shell_metadata <- function(x) {
    SummarizedExperiment(metadata=int_metadata(x))
}
