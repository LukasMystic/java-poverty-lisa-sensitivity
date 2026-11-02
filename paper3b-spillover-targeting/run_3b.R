###############################################################################
# Paper 3B - run the spillover-targeting pipeline.
#   Rscript run_3b.R        (from this directory)
#
# The scripts here resolve every path against their own location, so the two
# inputs built by the root script are copied in first if they are not present:
#   java_merged.gpkg          <- Rscript ../spatial_M1_M3.R
#   M3_lisaY_all_schemes.csv  <- same run
# Both copies are git-ignored; the originals stay at the repository root.
# Seed 20260826 throughout. Re-running overwrites every output in place, so
# `git diff` afterwards shows whether anything changed.
###############################################################################

here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else normalizePath(getwd())
})
setwd(here)
root <- dirname(here)

for (f in c("java_merged.gpkg", "M3_lisaY_all_schemes.csv")) {
  if (!file.exists(file.path(here, f))) {
    src <- file.path(root, f)
    if (!file.exists(src))
      stop("missing ", f, " - run `Rscript spatial_M1_M3.R` at the repository root first")
    file.copy(src, file.path(here, f))
    message("copied ", f, " from the repository root")
  }
}

steps <- c(
  "spatial_M4.R",      # first pass: spec tests, candidate models, SDM on AIC
  "spatial_M4b.R",     # re-examination under the written rule -> SDEM
  "spatial_M4c.R",     # theta robustness across eight graphs; rebuilds the figures
  "spatial_M5_M6.R"    # typology and the spillover-aware ranking
)

for (s in steps) {
  message("\n==== ", s, " ====")
  t <- system.time(source(file.path(here, s), local = new.env(), echo = FALSE))
  message(sprintf("   done in %.1f s", t["elapsed"]))
}
message("\nAll steps complete.")
