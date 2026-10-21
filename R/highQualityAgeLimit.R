## The oldest age that makeCohortDataFiles() keeps exact (its `minAgeForGrouping`); older ages are
## rounded to `pixelGroupAgeClass`, and -1 rounds every age.
##
## With `overrideAgeInFires`, ages up to `dataYear - earliestFireYear` are the ones the fire database
## can date (time since the last fire), so they are kept to the year. This does not depend on
## `overrideBiomassInFires`. It once did: with the biomass re-estimate off, rounding to 20-year
## classes turned ages 1-9 into 0, and Biomass_core never ages or grows an age-0 cohort.
## A `dataYear` before `earliestFireYear` leaves no fire-dated ages, so -1.
highQualityAgeLimit <- function(overrideAgeInFires, dataYear, earliestFireYear) {
  if (!isTRUE(overrideAgeInFires)) {
    return(-1)
  }
  limit <- dataYear - earliestFireYear
  if (isTRUE(limit >= 0)) limit else -1
}
