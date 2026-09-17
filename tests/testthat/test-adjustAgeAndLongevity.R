## With `adjustAgeAndLongevity = TRUE`, createBiomass_coreInputs() adjusted cohort ages right after
## building `pixelCohortData`, but the `LCCClassesToReplaceNN` block then rebuilt `pixelCohortData`
## from `pixelTable` and never adjusted it again: every run replacing a class (e.g. 240) carried
## unadjusted ages alongside the adjusted longevities in `sim$species`. This pins that each build
## of `pixelCohortData` is followed by an age adjustment before its biomass is partitioned. Static,
## so it needs no data.
test_that("every pixelCohortData build is age-adjusted before partitionBiomass", {
  code <- parse(file.path(moduleRoot, paste0(moduleName, ".R")))
  fn <- NULL
  for (expr in as.list(code)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("createBiomass_coreInputs"))) {
      fn <- expr[[3]]
    }
  }
  expect_false(is.null(fn))

  ## the calls of interest, in source order
  calls <- character(0)
  walk <- function(x) {
    if (is.call(x)) {
      nm <- if (is.name(x[[1]])) as.character(x[[1]]) else ""
      if (nm %in% c("makeAndCleanInitialCohortData", "adjustAgeToLongevity", "partitionBiomass")) {
        calls[[length(calls) + 1L]] <<- nm
      }
      lapply(as.list(x), walk)
    }
    invisible(NULL)
  }
  walk(fn)

  builds <- which(calls == "makeAndCleanInitialCohortData")
  expect_gte(length(builds), 1L)
  for (b in builds) {
    after <- calls[seq_along(calls) > b]
    nextPartition <- match("partitionBiomass", after)
    expect_false(is.na(nextPartition))
    expect_true("adjustAgeToLongevity" %in% after[seq_len(nextPartition - 1L)])
  }
})
