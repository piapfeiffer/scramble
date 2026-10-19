# Regenerates reference.rds, the reference values of the test
# "results match the reference values" in test-pcaSCRAMBLE.R.
#
# Only run this after a deliberate change of the numerical results, and
# document the change in NEWS.md. Run from the package root:
#   Rscript tests/testthat/fixtures/make-reference.R
devtools::load_all(quiet = TRUE)
source("tests/testthat/helper.R")
x <- make_data()
ref <- lapply(reference_configs, function(args) {
  do.call(fit_quick, c(list(x), args))[reference_fields]
})
saveRDS(ref, "tests/testthat/fixtures/reference.rds", version = 2)
