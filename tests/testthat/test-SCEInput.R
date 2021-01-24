# This tests the SCEInput class and related methods.
# library(testthat); library(SingleCellExperiment); source("setup.R"); source("test-SCEInput.R")

test_that("MainExpInput works as expected", {
    out <- MainExpInput(a=1, b="2")
    expect_identical(inputArguments(out), list(a=1, b="2"))
    expect_identical(loaded, getInput(loaded, out))

    expect_error(show(out), NA)

    expect_true(isInputSameShape(out))
    expect_false(isInputReducedDim(out))
    expect_false(isInputAltExp(out))
})

test_that("AssayInput works as expected", {
    out <- AssayInput(assay="logcounts", a=1, b="2")
    expect_identical(inputArguments(out), list(a=1, b="2"))
    expect_identical(logcounts(loaded), getInput(loaded, out))

    expect_error(show(out), NA)
    expect_error(AssayInput(assay=letters), "scalar")
    expect_error(AssayInput(assay=0i), "string or integer")

    expect_true(isInputSameShape(out))
    expect_false(isInputReducedDim(out))
    expect_false(isInputAltExp(out))
})

test_that("ReducedDimInput works as expected", {
    out <- ReducedDimInput(type="TSNE", a=1, b="2")
    expect_identical(inputArguments(out), list(a=1, b="2"))
    expect_identical(reducedDim(loaded, "TSNE"), getInput(loaded, out))

    expect_error(show(out), NA)
    expect_error(ReducedDimInput(type=letters), "scalar")
    expect_error(ReducedDimInput(type=0i), "string or integer")

    expect_false(isInputSameShape(out))
    expect_true(isInputReducedDim(out))
    expect_false(isInputAltExp(out))
})

test_that("AltExpInput works as expected", {
    out <- AltExpInput(experiment="Protein", a=1, b="2")
    expect_identical(inputArguments(out), list(a=1, b="2"))
    expect_identical(altExp(loaded, 2), getInput(loaded, out))

    expect_error(show(out), NA)
    expect_error(AltExpInput(experiment=0i), "string or integer")

    expect_true(isInputSameShape(out))
    expect_false(isInputReducedDim(out))
    expect_true(isInputAltExp(out))
})

test_that("AltAssayInput works as expected", {
    out <- AltAssayInput(experiment="Spike", assay="logcounts", a=1, b="2")
    expect_identical(inputArguments(out), list(a=1, b="2"))
    expect_identical(assay(altExp(loaded), "logcounts"), getInput(loaded, out))

    expect_error(show(out), NA)
    expect_error(AltExpInput(experiment=0i), "string or integer")
    expect_error(AltAssayInput(experiment=letters, assay=2), "scalar")
    expect_error(AltAssayInput(experiment=2, assay=letters), "scalar")

    expect_true(isInputSameShape(out))
    expect_false(isInputReducedDim(out))
    expect_true(isInputAltExp(out))
})

test_that("AltReducedDimInput works as expected", {
    altExp(loaded, 2) <- as(altExp(loaded, 2), "SingleCellExperiment")
    reducedDim(altExp(loaded, 2), "TSNE") <- reducedDim(loaded) * 10

    out <- AltReducedDimInput(experiment=2, type="TSNE", a=1, b="2")
    expect_identical(inputArguments(out), list(a=1, b="2"))
    expect_identical(reducedDim(altExp(loaded, 2), "TSNE"), getInput(loaded, out))

    expect_error(show(out), NA)
    expect_error(AltExpInput(experiment=0i), "string or integer")
    expect_error(AltReducedDimInput(experiment=2, type=letters), "scalar")
    expect_error(AltReducedDimInput(experiment=2, type=0i), "string or integer")

    expect_false(isInputSameShape(out))
    expect_true(isInputReducedDim(out))
    expect_true(isInputAltExp(out))
})
