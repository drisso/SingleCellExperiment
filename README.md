---
title: "An introduction to the SingleCellExperiment class"
author: "Davide Risso"
date: "Last modified: June 6, 2017; Compiled: `r format(Sys.time(), '%B %d, %Y')`"
output:
  BiocStyle::html_document:
    toc: true
vignette: >
  %\VignetteEncoding{UTF-8}
---

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{SingleCellExperiment Vignette}
-->

```{r options, include=FALSE, echo=FALSE}
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)
```

# Package Status

|                |               |
| -------------- | ------------- |
| Project Status | [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) |
| Travis CI      | [![Build Status](https://travis-ci.org/drisso/SingleCellExperiment.svg?branch=master)](https://travis-ci.org/drisso/SingleCellExperiment) |
| Test coverage  | [![Coverage](https://codecov.io/gh/drisso/SingleCellExperiment/branch/master/graph/badge.svg)](https://codecov.io/gh/drisso/SingleCellExperiment) |

# Setup

This vignette requires the [SingleCellExperiment](https://github.com/drisso/SingleCellExperiment) package,
available from GitHub, and the `scRNAseq` package, available from Bioconductor.

```{r install, eval=FALSE}
BiocInstaller::biocLite("drisso/SingleCellExperiment")
BiocInstaller::biocLite("scRNAseq")
```

```{r load_packages}
library(SingleCellExperiment)
library(scRNAseq)
```

# The SingleCellExperiment class

The `SingleCellExperiment` class is a light-weight container for single-cell 
genomics data that extends the `RangedSummarizedExperiment` class with the following
additional slots and methods specific to single-cell genomics datasets.

* `int_elementMetadata` 
* `int_colData`
* `int_metadata`
* `reducedDims`

As suggested by the `int_` prefix, the first three slots are not meant for direct
manipulation, but rather are set by other methods. For instance, `isSpike<-` will
set a proper column of `int_elementMetadata` and `sizeFactors<-` will set a 
column in the `int_colData` slot.

## Create instances of SingleCellExperiment

There are two main ways to create instances of `SingleCellExperiment`. The first
is via the constructor.

```{r construct}
sce <- SingleCellExperiment(
  assays = list(counts = matrix(rpois(100, lambda = 10), ncol=10, nrow=10)))
sce
```

The second is via coercion from `SummarizedExperiment` objects.

```{r coerce}
se <- SummarizedExperiment(
  assays = list(counts = matrix(rpois(100, lambda = 10), ncol=10, nrow=10)))
as(se, "SingleCellExperiment")
```

# A simple example

Here we use a subset of the `allen` dataset from the `scRNAseq` package to
show a typical use of the class.

The `allen` data are stored as a `SummarizedExperiment`, so we will use the
coercion method to turn it into a `SingleCellExperiment`.

```{r fluidigm}
data(allen)
allen

sce <- as(allen, "SingleCellExperiment")
sce
```

## Adding spike-in information

One of the main additions to `SummarizedExperiment` is the ability for the user
to add information about the spike-ins, with the method `isSpike`.

```{r spikes}
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
sce
table(isSpike(sce))
spikeNames(sce)
```

Although for the majority of cases one set of spike-ins should be enough, the
class has the flexibility of including more than one set of spikes.
Let us pretend that the members of the Adam gene family have beeen spiked-in as
external genes in these data.

```{r spikes2}
isSpike(sce, "Adam") <- grepl("^Adam[0-9]", rownames(sce))
sce
table(isSpike(sce))
table(isSpike(sce, "ERCC"))
table(isSpike(sce, "Adam"))
spikeNames(sce)
```

## Adding size factors

Similarly, one can include the computed size factors. For illustration, we simply
compute the total number of reads as size factors, but better ways to compute 
size factors are available (see, e.g., the [scran](https://www.bioconductor.org/packages/scran) package).

```{r sizeFactors}
sizeFactors(sce) <- colSums(assay(sce))
head(sizeFactors(sce))
```

We can compute multiple size factors and store them in the object, by providing
a name.

```{r sizeFactors2}
sizeFactors(sce, "ERCC") <- colSums(assay(sce)[isSpike(sce, "ERCC"),])
head(sizeFactors(sce, "ERCC"))
```

## Retrieve `colData` and `rowData` information

The `colData` and `rowData` methods can be used to retrieve the stored sample-
and gene-level metadata. By default, spike-ins and size factors are not returned
by such methods, as they are conceptually distinct from the rest of the metadata.

```{r metadata}
colData(sce)
rowData(sce)
```

However, it is sometimes useful to retrieve a `DataFrame` with all the available
metadata. This can be achieved by specifying `internal=TRUE`.

```{r metadata2}
colData(sce, internal=TRUE)
rowData(sce, internal=TRUE)
```

## Adding low-dimentional representations

For simplicity and speed we work on a subset of 100 genes. To avoid
ending up with only uninteresting genes, we extract the 100 genes with maximal 
variance.

```{r subset}
library(magrittr)
assay(sce) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)

sce_sub <- sce[names(vars[1:100]),]
sce_sub
```

We then obtain the PCA and t-SNE representations of the data and add them to the
object with the `reducedDims` method

```{r pca}
library(Rtsne)
set.seed(5252)

pca_data <- prcomp(t(log1p(assay(sce_sub))))
tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)

reducedDims(sce_sub) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data$Y)
sce_sub

reducedDims(sce_sub)
head(reducedDim(sce_sub, "PCA")[,1:2])
head(reducedDim(sce_sub, "TSNE")[,1:2])
```

# Design decisions

In principle, all we want from the additional slot is a mechanism to protect
some row and column data from the user direct manipulation, such that, for
instance, only a call to `sizeFactors<-` can change the size factors.
If there was a way to reserve a subset of columns (or column names) as "private"
in `colData()` and `rowData()` we would not need the additional slots.

For the `reducedDims` slot, things are slightly different, since one can imagine
that different dimentionality techniques will be useful for different aspects
of the analysis (e.g., t-SNE for visualization, PCA for pseudo-time inference,
etc.). We see `reducedDims` as a similar slot to `assays()`.

We decided to extend `RangedSummarizedExperiment` rather than `SummarizedExperiment`
because for certain assays it will be essential to have `rowRanges()`. Even for
RNA-seq it is sometimes useful to have `rowRanges()` and other classes, e.g.,
`DESeqDataSet` in the `DESeq2` package.
An alternative would have been to have two classes, `SingleCellExperiment` and
`RangedSingleCellExperiment`, but this seems an unnecessary duplication as
having a class with default empty rowRanges seems good enough when one does not
need rowRanges.

By design, the scope of this package is rather limited: i.e., to define the
`SingleCellExperiment` class and some minimal getter and setter methods.
For this reason, we leave it to developers to provide more advance methods for
the `SingleCellExperiment` class in specialized packages that depend on this.
For instance, it would be natural to write a `pca` and a `tsne` method for the 
class. Similarly, we leave it to the developers of the existing single-cell
packages to provide coercion methods from their specialized class to 
`SingleCellExperiment`.

# Session Info

```{r}
sessionInfo()
```
