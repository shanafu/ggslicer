# Locate the real-image testdata/ fixture directory (a sibling of every repo,
# never committed to git). Honors GGSLICER_TESTDATA_DIR if set; otherwise
# tries a few relative candidates, since testthat's working directory during
# a test run can be either the package root or tests/testthat depending on
# how tests are invoked (devtools::test() vs. R CMD check vs. testthat::test_file()).
testdata_dir <- function() {
  env_dir <- Sys.getenv("GGSLICER_TESTDATA_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(env_dir)
  }
  candidates <- c("../testdata", "../../testdata", "../../../testdata")
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[1])
  }
  candidates[1]
}

# Skip the current test if testdata/ isn't available on this machine.
skip_if_no_testdata <- function() {
  testthat::skip_if_not(dir.exists(testdata_dir()), "testdata/ not available")
}
