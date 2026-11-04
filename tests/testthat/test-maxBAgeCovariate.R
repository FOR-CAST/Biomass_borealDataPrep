## `standAgeMapForMaxB` replaces the age covariate of the maxB fit, and only that: the cohorts, their
## biomass and cover, the strata and cohortData keep `standAgeMap`.
library(data.table)

moduleBody <- function(name) {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(name)), exprs)
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

callsTo <- function(expr, fn) {
  out <- list()
  if (!is.call(expr)) return(out)
  if (identical(expr[[1]], as.name(fn))) out <- list(expr)
  for (i in seq_along(expr)) {
    if (is.call(expr[[i]])) out <- c(out, callsTo(expr[[i]], fn))
  }
  out
}

## 2 x 3 grid; cell 2 has no age, cell 3 age 0, and cell 4's 0.4 rounds to 0
ageMap <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2,
                      crs = "EPSG:3978", vals = c(80, NA, 0, 0.4, 120, 45))
## the cohorts' own age (standAgeMap), cell by cell; makePixelTable() logs cell 2's 30.4 unrounded
standAgeMap <- terra::setValues(ageMap, c(100, 30.4, 5, 60, 140, 45))
cohorts <- data.table(
  pixelIndex = c(1L, 1L, 2L, 3L, 4L, 5L, 6L),
  speciesCode = factor(c("Pice_mar", "Popu_tre", "Pice_mar", "Pice_mar", "Pinu_ban", "Pice_gla", "Pice_mar")),
  ecoregionGroup = factor(c("1_210", "1_210", "1_220", "1_220", "2_210", "2_210", "2_810")),
  age = c(100L, 100L, 30L, 5L, 60L, 140L, 45L),
  logAge = log(c(100, 100, 30, 5, 60, 140, 45)),
  B = c(5000, 2000, 3000, 1000, 4000, 9000, 2500),
  cover = c(60L, 40L, 100L, 100L, 100L, 100L, 100L),
  lcc = c(210L, 210L, 220L, 220L, 210L, 210L, 810L),
  totalBiomass = c(7000L, 7000L, 3000L, 1000L, 4000L, 9000L, 2500L)
)
default <- unique(cohorts[, .(B, logAge, speciesCode, ecoregionGroup, lcc, cover)])

test_that("maxBModelData() takes logAge from the map and leaves everything else alone", {
  ## 1_220 has no age in the map in any cohort: see the fallback tests below
  out <- maxBModelData(cohorts, ageMap, standAgeMap)[ecoregionGroup != "1_220"]
  expect_identical(names(out), names(default))
  kept <- cohorts$pixelIndex %in% c(1L, 5L, 6L)
  expect_identical(out$B, cohorts$B[kept])
  expect_identical(out$cover, cohorts$cover[kept])
  expect_identical(out$speciesCode, cohorts$speciesCode[kept])
  expect_identical(out$ecoregionGroup, cohorts$ecoregionGroup[kept])
  expect_identical(out$lcc, cohorts$lcc[kept])
  expect_equal(out$logAge, log(c(80, 80, 120, 45)))
})

test_that("maxBModelData() leaves out a cohort with no age in the map where its group has one", {
  ## each group has an aged cohort: no age (cell 2) and age 0 (cell 4) are left out, though both
  ## cohorts have an age in standAgeMap; so is cell 4's 0.4, which rounds to 0, in the first map
  am <- terra::setValues(ageMap, c(80, NA, 25, 0, 120, 45))
  out <- maxBModelData(cohorts, am, standAgeMap)
  expect_identical(out$B, cohorts$B[cohorts$pixelIndex %in% c(1L, 3L, 5L, 6L)])
  expect_equal(out$logAge, log(c(80, 80, 25, 120, 45)))
  expect_identical(attr(out, "nDropped"), 2L)
  expect_length(attr(out, "nFallback"), 0L)

  fit <- maxBModelData(cohorts, ageMap, standAgeMap)
  expect_false(any(fit$speciesCode == "Pinu_ban"))
  expect_identical(attr(fit, "nDropped"), 1L)
})

test_that("maxBModelData() fits a group with no age in the map at its cohorts' age in standAgeMap", {
  out <- maxBModelData(cohorts, ageMap, standAgeMap)
  fb <- out[ecoregionGroup == "1_220"]
  expect_identical(fb$B, c(3000, 1000))
  expect_equal(fb$logAge, log(c(30.4, 5)))
  expect_identical(as.character(fb$speciesCode), c("Pice_mar", "Pice_mar"))
  expect_identical(fb$lcc, c(220L, 220L))
  expect_identical(fb$cover, c(100L, 100L))
  expect_identical(names(out), names(default))
  expect_identical(attr(out, "nFallback"), c(`1_220` = 2L))

  ## the same rule there: no age, or age 0, in standAgeMap is left out
  sam <- terra::setValues(standAgeMap, c(100, 0, 5, 60, 140, 45))
  out0 <- maxBModelData(cohorts, ageMap, sam)
  expect_equal(out0[ecoregionGroup == "1_220", logAge], log(5))
  expect_identical(attr(out0, "nFallback"), c(`1_220` = 1L))
  expect_identical(attr(out0, "nDropped"), 2L)
})

test_that("a group with no age in either map is not fitted, and its speciesEcoregion rows stop", {
  sam <- terra::setValues(standAgeMap, c(100, NA, 0, 60, 140, 45))
  fit <- maxBModelData(cohorts, ageMap, sam)
  expect_false(any(fit$ecoregionGroup == "1_220"))
  expect_identical(attr(fit, "nFallback"), c(`1_220` = 0L))
  expect_identical(attr(fit, "nDropped"), 3L)
  expect_error(maxBSpeciesEcoregionRows(fit, default),
               "neither standAgeMapForMaxB nor standAgeMap has an age .* ecoregionGroup\\(s\\) 1_220")
})

test_that("maxBModelData() floors logAge as makePixelTable() does and keeps the default's unique()", {
  am <- terra::setValues(ageMap, c(0.6, 50, 50, 50, 50, 50))
  twins <- data.table::copy(cohorts)[c(5L, 6L)]
  twins[, `:=`(pixelIndex = c(2L, 3L), speciesCode = factor("Pinu_ban"), ecoregionGroup = factor("2_210"),
               B = 4000, cover = 100L, lcc = 210L)]
  out <- maxBModelData(rbind(cohorts[1L], twins), am, standAgeMap)
  expect_equal(out$logAge, c(log(0.6), log(50)))
  expect_identical(nrow(out), 2L)
})

test_that("maxBModelData() refuses a map on a smaller grid than the cohorts' pixelIndex", {
  small <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                       crs = "EPSG:3978", vals = 50)
  expect_error(maxBModelData(cohorts, small, standAgeMap), "rasterToMatch_biomassParam grid")
  ## and so is standAgeMap, read for 1_220 (cells 2 and 3)
  one <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                     crs = "EPSG:3978", vals = 50)
  expect_error(maxBModelData(cohorts, ageMap, one), "^standAgeMap has 1 cells")
})

test_that("maxBSpeciesEcoregionRows() keeps the default fit's species x ecoregionGroup rows", {
  ## no cohort of 1_220 has an age in either map
  fit <- maxBModelData(cohorts, ageMap, terra::setValues(standAgeMap, c(100, NA, 0, 60, 140, 45)))
  ## nothing missing: the rows come back unchanged
  same <- maxBSpeciesEcoregionRows(default, default)
  expect_identical(same, default)
  expect_identical(attr(same, "nAdded"), 0L)

  ## so the fit knows nothing of that group: stop
  expect_error(maxBSpeciesEcoregionRows(fit, default), "ecoregionGroup\\(s\\) 1_220")

  ## with 1_220 in the fit through another species, its Pice_mar row is added back, and so is
  ## Pinu_ban in 2_210 (the fit has 2_210 through Pice_gla)
  fit2 <- rbind(fit, data.table(B = 1500, logAge = log(20), speciesCode = factor("Popu_tre"),
                                ecoregionGroup = factor("1_220"), lcc = 220L, cover = 100L))
  rows <- maxBSpeciesEcoregionRows(fit2, default)
  expect_identical(attr(rows, "nAdded"), 2L)
  key <- c("ecoregionGroup", "speciesCode")
  expect_true(all(unique(default[, ..key])[, paste(ecoregionGroup, speciesCode)] %in%
                    rows[, paste(ecoregionGroup, speciesCode)]))
  expect_identical(names(rows), names(fit2))
})

test_that("in createBiomass_coreInputs() only the maxB fit reads standAgeMapForMaxB", {
  bod <- moduleBody("createBiomass_coreInputs")
  ## two statements name it: the alignment to rasterToMatch_biomassParam and the fit's rows
  mentions <- Filter(function(x) "standAgeMapForMaxB" %in% all.names(x), as.list(bod))
  expect_length(mentions, 2L)
  expect_true(any(vapply(mentions, function(x) length(callsTo(x, "postProcess")) > 0, logical(1))))
  expect_true(any(vapply(mentions, function(x) length(callsTo(x, "maxBModelData")) > 0, logical(1))))
  ## a group with no age in it falls back to standAgeMap, the cohorts' own age
  mbd <- callsTo(bod, "maxBModelData")
  expect_length(mbd, 1L)
  expect_identical(mbd[[1]]$fallbackAgeMap, quote(sim$standAgeMap))

  ## the cohorts are built from standAgeMap
  mpt <- callsTo(bod, "makePixelTable")
  expect_length(mpt, 1L)
  expect_identical(mpt[[1]]$standAgeMap, quote(sim$standAgeMap))

  ## the biomass model's subsample is drawn from the fit's rows ...
  sdt <- Filter(function(x) identical(x[[2]], as.name("biomassModelData")), callsTo(bod, "subsetDT"))
  expect_length(sdt, 1L)
  ## ... and the observed-biomass floor still uses the default rows
  fm <- callsTo(bod, "floorMaxBAtObserved")
  expect_length(fm, 1L)
  expect_identical(fm[[1]][[3]], as.name("cohortDataOnlyForestLCCBiomass"))
})

test_that("without standAgeMapForMaxB the fit keeps its default rows, cache key and speciesEcoregion rows", {
  bod <- moduleBody("createBiomass_coreInputs")

  ## NULL: the fit's rows are the default ones and nothing is added to its cache key
  br <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                 identical(x[[2]], quote(is.null(sim$standAgeMapForMaxB))), as.list(bod))
  expect_length(br, 1L)
  expect_identical(as.list(br[[1]][[3]])[-1L],
                   list(quote(biomassModelData <- cohortDataOnlyForestLCCBiomass),
                        quote(maxBAgeKey <- NULL)))

  ## each cache key of the fit ends with maxBAgeKey, which c() drops when it is NULL
  keys <- Filter(function(x) identical(x[[2]], as.name("sumResponse")), callsTo(bod, "<-"))
  expect_length(keys, 2L)
  for (k in keys) {
    rhs <- k[[3]]
    expect_identical(rhs[[1]], as.name("c"))
    expect_identical(rhs[[length(rhs)]], as.name("maxBAgeKey"))
  }
  expect_identical(c(1234, NULL), 1234)

  ## speciesEcoregion is made from the fit's subsample, as before the input existed
  ser <- Filter(function(x) identical(x[[2]], as.name("speciesEcoregionRows")), callsTo(bod, "<-"))
  expect_length(ser, 1L)
  rhs <- ser[[1]][[3]]
  expect_identical(rhs[[2]], quote(is.null(maxBAgeKey)))
  expect_identical(rhs[[3]][[1]], as.name("{"))
  expect_identical(as.list(rhs[[3]])[-1L], list(as.name("cohortDataOnlyForestLCCBiomassSubset")))
})
