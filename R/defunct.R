#' Defunct methods
#'
#' Defunct methods in the \pkg{SingleCellExperiment} package.
#'
#' @section Named size factors:
#' The class now only supports one set of size factors, accessible via \code{\link{sizeFactors}}.
#' This represents a simplification of the class and removes a difficult part of the API 
#' (that had to deal with both \code{NULL} and strings to specify the size factor set of interest).
#'
#' @section Spike-ins:
#' It is recommended to handle spike-ins and other \dQuote{alternative} features via \code{\link{altExps}}.
#'
#' @author Aaron Lun
#'
#' @docType methods 
#' @aliases isSpike isSpike<- spikeNames clearSpikes
#' sizeFactorNames clearSizeFactors
#' @name defunct
NULL
