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

.filled_int_colData <- function(x) {
    df <- int_colData(x)
    .fill_int_fields(df, c(.red_key, .alt_key))
}
