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
  out <- maxBModelData(cohorts, ageMap)
  expect_identical(names(out), names(default))
  kept <- cohorts$pixelIndex %in% c(1L, 5L, 6L)
  expect_identical(out$B, cohorts$B[kept])
  expect_identical(out$cover, cohorts$cover[kept])
  expect_identical(out$speciesCode, cohorts$speciesCode[kept])
  expect_identical(out$ecoregionGroup, cohorts$ecoregionGroup[kept])
  expect_identical(out$lcc, cohorts$lcc[kept])
  expect_equal(out$logAge, log(c(80, 80, 120, 45)))
  ## no age (cell 2), age 0 (cell 3), and an age that rounds to 0 (cell 4)
  expect_identical(attr(out, "nDropped"), 3L)
})

test_that("maxBModelData() floors logAge as makePixelTable() does and keeps the default's unique()", {
  am <- terra::setValues(ageMap, c(0.6, 50, 50, 50, 50, 50))
  twins <- data.table::copy(cohorts)[c(5L, 6L)]
  twins[, `:=`(pixelIndex = c(2L, 3L), speciesCode = factor("Pinu_ban"), ecoregionGroup = factor("2_210"),
               B = 4000, cover = 100L, lcc = 210L)]
  out <- maxBModelData(rbind(cohorts[1L], twins), am)
  expect_equal(out$logAge, c(log(0.6), log(50)))
  expect_identical(nrow(out), 2L)
})

test_that("maxBModelData() refuses a map on a smaller grid than the cohorts' pixelIndex", {
  small <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                       crs = "EPSG:3978", vals = 50)
  expect_error(maxBModelData(cohorts, small), "rasterToMatch_biomassParam grid")
})

test_that("maxBSpeciesEcoregionRows() keeps the default fit's species x ecoregionGroup rows", {
  fit <- maxBModelData(cohorts, ageMap)
  ## nothing missing: the rows come back unchanged
  same <- maxBSpeciesEcoregionRows(default, default)
  expect_identical(same, default)
  expect_identical(attr(same, "nAdded"), 0L)

  ## no cohort of 1_220 has an age in the map, so the fit knows nothing of that group: stop
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
