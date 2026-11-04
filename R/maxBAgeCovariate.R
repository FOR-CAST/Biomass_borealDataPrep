## The age covariate of `biomassModel`, the fit that gives maxB and maxANPP.
##
## By default the fit takes each cohort's age from `standAgeMap`, the age the cohorts start the
## simulation with. A project can supply `standAgeMapForMaxB` to fit maxB and maxANPP with another
## age layer -- e.g. the age product made with the biomass map the fit's B comes from -- while the
## cohorts keep `standAgeMap`. Only the fit's `logAge` changes: the cohorts, their biomass and cover,
## the strata, `cohortData` and the prediction at `cover = 100` and age = longevity are the same.
## An ecoregionGroup with no age in `standAgeMapForMaxB` in any cohort is fitted with
## `standAgeMap`'s (see maxBModelData()).

## The biomass model's rows, with `logAge` from `ageMap` at each cohort's pixel.
##
## `cohortData` is the module's estimation cohorts with biomass (`totalBiomass > 0`), `pixelIndex`
## on the grid of `ageMap` (`rasterToMatch_biomassParam`). A cohort whose pixel has no age in
## `ageMap`, or age 0, is left out: the module gives an age-0 stand no biomass
## (`LandR::makeAndCleanInitialCohortData()`), so its biomass and age would disagree. An
## ecoregionGroup in which no cohort has an age in `ageMap` would have no rows, and maxB could not
## be predicted there (`LandR::makeSpeciesEcoregion()` knows only the fit's groups), so its cohorts
## take their age from `fallbackAgeMap` (`standAgeMap`, the age the default path fits them with),
## under the same rule. `logAge` is computed as `LandR::makePixelTable()` computes it. Returns the
## unique rows, as the default path does, with attributes "nDropped": the cohort rows left out, and
## "nFallback": the cohort rows each fallback ecoregionGroup gives the fit, named by group.
maxBModelData <- function(cohortData, ageMap, fallbackAgeMap) {
  ageForFit <- maxBCellAge(ageMap, cohortData$pixelIndex, "standAgeMapForMaxB")
  group <- as.character(cohortData$ecoregionGroup)
  noAgeGroups <- sort(setdiff(group, group[maxBHasAge(ageForFit)]))
  fallback <- group %in% noAgeGroups
  if (any(fallback)) {
    ageForFit[fallback] <- maxBCellAge(fallbackAgeMap, cohortData$pixelIndex[fallback], "standAgeMap")
  }
  hasAge <- maxBHasAge(ageForFit)
  nFallback <- vapply(noAgeGroups, function(g) sum(hasAge[group == g]), integer(1))
  keepRows <- which(hasAge)
  ## the columns, and their order, of the default path's `.(B, logAge, speciesCode, ...)`; built
  ## from vectors because `cohortData` has its own `age` column
  out <- data.table::data.table(
    B = cohortData$B[keepRows],
    logAge = LandR:::.logFloor(ageForFit[keepRows]),
    speciesCode = cohortData$speciesCode[keepRows],
    ecoregionGroup = cohortData$ecoregionGroup[keepRows],
    lcc = cohortData$lcc[keepRows],
    cover = cohortData$cover[keepRows]
  )
  out <- unique(out)
  data.table::setattr(out, "nDropped", nrow(cohortData) - length(keepRows))
  data.table::setattr(out, "nFallback", nFallback)
  out
}

## The age in `ageMap` at cells `pixelIndex`; `name` is the input it is, for the error.
maxBCellAge <- function(ageMap, pixelIndex, name) {
  vals <- terra::values(ageMap, mat = FALSE)
  if (any(pixelIndex > length(vals))) {
    stop(name, " has ", length(vals), " cells, fewer than the cohorts' pixelIndex; ",
         "it must be on the rasterToMatch_biomassParam grid")
  }
  vals[pixelIndex]
}

## An age the fit can use: present, and not 0 once rounded as cohort ages are.
maxBHasAge <- function(age) {
  !is.na(age) & LandR::asInteger(age) > 0L
}

## The rows `LandR::makeSpeciesEcoregion()` makes speciesEcoregion from.
##
## It takes one row per species x ecoregionGroup present in the data it is given -- the model's
## subsample -- and predicts maxB for each. When the fit uses `standAgeMapForMaxB`, a species x
## ecoregionGroup whose cohorts all lack an age there would lose its row, and then its cohorts
## (makeCohortDataFiles() keeps only those with parameters). One row of each such combination from
## `allRows`, the default fit's data, is added, so speciesEcoregion has the rows it has without
## `standAgeMapForMaxB`; its maxB is predicted from the fit's species and ecoregionGroup effects.
## An ecoregionGroup the fit has no data for -- no cohort there has an age in either age map (see
## maxBModelData()) -- cannot be predicted, and stops.
## Returns `fitRows` unchanged when nothing is missing; attribute "nAdded" counts the added rows.
maxBSpeciesEcoregionRows <- function(fitRows, allRows) {
  key <- c("ecoregionGroup", "speciesCode")
  have <- unique(fitRows[, key, with = FALSE])
  missingRows <- unique(allRows, by = key)[!have, on = key]
  if (nrow(missingRows) == 0L) {
    data.table::setattr(fitRows, "nAdded", 0L)
    return(fitRows)
  }
  noData <- setdiff(as.character(missingRows$ecoregionGroup), as.character(fitRows$ecoregionGroup))
  if (length(noData)) {
    stop("neither standAgeMapForMaxB nor standAgeMap has an age in any cohort of ecoregionGroup(s) ",
         paste(sort(unique(noData)), collapse = ", "), ", so maxB cannot be fitted there")
  }
  out <- data.table::rbindlist(list(fitRows, missingRows[, names(fitRows), with = FALSE]))
  data.table::setattr(out, "nAdded", nrow(missingRows))
  out
}
