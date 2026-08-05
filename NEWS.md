Known issues: <https://github.com/PredictiveEcology/Biomass_borealDataPrep/issues>

version 1.5.13
=============

## dependency changes
* requires `LandR (>= 1.2.0.9005)`, which reworked `convertUnwantedLCC()`: it gains a
  `method` argument, and it again returns the `newPossLCC` column this module needs to
  write replacement classes back into `rstLCCAdj`. That column was absent in LandR
  1.2.0.9004, where the `is.null()` guard below silently stopped firing and left
  `rstLCCAdj` (hence `ecoregionMap`) showing the un-replaced classes.

## new features
* new parameter `LCCClassesToReplaceNNMethod`, passed to `LandR::convertUnwantedLCC()` as
  its `method`. Default `"nearestWeighted"` allocates each `LCCClassesToReplaceNN` pixel a
  neighbouring class drawn in proportion to that class's local abundance, keyed on the
  pixel's ground position so the result is deterministic without a seed and a grid-aligned
  crop of the study area agrees with the full extent. `"nearestRandom"` draws from the RNG
  instead, for when replicates should differ.

version 1.5.4
=============

## dependency changes
* terra is now a dependency

## new features


## bug fixes
* Comment fixes
