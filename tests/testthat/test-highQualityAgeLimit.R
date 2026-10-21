## Ages the fire database can date are kept exact by makeCohortDataFiles(), whatever
## overrideBiomassInFires is. The limit used to be set only inside the overrideBiomassInFires block,
## so with the biomass re-estimate off every age was rounded to pixelGroupAgeClass, and ages 1-9
## became 0, an age Biomass_core does not advance until its first age reclassification.

moduleBody <- function(name) {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(name)), exprs)
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

## every call to `fn` anywhere in `expr` (a `|>` pipe is already a nested call once parsed)
callsTo <- function(expr, fn) {
  out <- list()
  if (!is.call(expr)) return(out)
  if (identical(expr[[1]], as.name(fn))) out <- list(expr)
  for (i in seq_along(expr)) {
    if (is.call(expr[[i]])) out <- c(out, callsTo(expr[[i]], fn))
  }
  out
}

## every `name <- ...` anywhere in `expr`
assignmentsTo <- function(expr, name) {
  Filter(function(x) identical(x[[2]], as.name(name)), callsTo(expr, "<-"))
}

test_that("highQualityAgeLimit() is dataYear - earliestFireYear with overrideAgeInFires, else -1", {
  expect_identical(highQualityAgeLimit(TRUE, 2020, 1950L), 70)
  expect_identical(highQualityAgeLimit(FALSE, 2020, 1950L), -1)
  expect_identical(highQualityAgeLimit(NULL, 2020, 1950L), -1)
  expect_identical(highQualityAgeLimit(TRUE, 1940, 1950L), -1) ## no fire-dated ages
  expect_identical(highQualityAgeLimit(TRUE, 1950, 1950L), 0)
})

test_that("createBiomass_coreInputs() sets the limit outside the overrideBiomassInFires block", {
  bod <- moduleBody("createBiomass_coreInputs")

  ## one assignment, a top-level statement of the function, from highQualityAgeLimit()
  expect_length(assignmentsTo(bod, "maxAgeHighQualityData"), 1L)
  top <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("maxAgeHighQualityData")), as.list(bod))
  expect_length(top, 1L)
  rhs <- top[[1]][[3]]
  expect_identical(rhs[[1]], as.name("highQualityAgeLimit"))
  expect_identical(rhs$overrideAgeInFires, quote(P(sim)$overrideAgeInFires))
  expect_identical(rhs$dataYear, quote(P(sim)$dataYear))
  expect_identical(rhs$earliestFireYear, quote(P(sim)$earliestFireYear))

  ## the overrideBiomassInFires block reads it and no longer sets it
  overrideBlock <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                            identical(x[[2]], quote(isTRUE(P(sim)$overrideBiomassInFires))),
                          as.list(bod))
  expect_length(overrideBlock, 1L)
  expect_length(assignmentsTo(overrideBlock[[1]], "maxAgeHighQualityData"), 0L)

  ## and it is what makeCohortDataFiles() gets as minAgeForGrouping
  mcdf <- callsTo(bod, "makeCohortDataFiles")
  expect_length(mcdf, 1L)
  expect_identical(mcdf[[1]]$minAgeForGrouping, as.name("maxAgeHighQualityData"))
})

test_that("with the limit makeCohortDataFiles() keeps ages up to it exact; without it 1-9 become 0", {
  withr::local_options(reproducible.useCache = FALSE)
  ages <- c(3L, 9L, 10L, 35L, 70L, 95L)
  pcd <- data.table::data.table(
    pixelIndex = seq_along(ages), speciesCode = factor("Pice_mar"), ecoregionGroup = factor("1_210"),
    age = ages, B = 1000L * seq_along(ages), totalBiomass = 1000L * seq_along(ages), cover = 100L,
    lcc = 210L, logAge = log(ages)
  )
  se <- data.table::data.table(ecoregionGroup = factor("1_210"), speciesCode = factor("Pice_mar"),
                               maxB = 10000L, maxANPP = 333L, establishprob = 0.5)
  run <- function(limit) {
    out <- suppressMessages(LandR::makeCohortDataFiles(
      pixelCohortData = data.table::copy(pcd),
      columnsForPixelGroups = LandR::columnsForPixelGroups(),
      speciesEcoregion = se,
      pixelGroupBiomassClass = 100,
      pixelGroupAgeClass = 20,
      minAgeForGrouping = limit,
      rmImputedPix = FALSE,
      imputedPixID = integer(0),
      pixelFateDT = LandR::pixelFate(fate = "Total number pixels", runningPixelTotal = length(ages))
    ))
    out$pixelCohortData[order(pixelIndex), age]
  }

  ## overrideAgeInFires = TRUE, dataYear 2020, earliestFireYear 1950: ages to 70 kept, 95 -> 100
  expect_identical(run(highQualityAgeLimit(TRUE, 2020, 1950L)), c(3L, 9L, 10L, 35L, 70L, 100L))
  ## every age rounded to 20 years
  expect_identical(run(highQualityAgeLimit(FALSE, 2020, 1950L)), c(0L, 0L, 20L, 40L, 80L, 100L))
})
