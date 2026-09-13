test_that("check_sitk_image() passes a real SimpleITK image through silently", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(2, 2, 2), function(idx) 0)
  expect_true(isTRUE(check_sitk_image(img)))
})

test_that("check_sitk_image() gives a specific, actionable error for a file path", {
  err <- tryCatch(check_sitk_image("brain.nii.gz"), error = function(e) conditionMessage(e))
  expect_match(err, "not a file path")
  expect_match(err, "brain.nii.gz", fixed = TRUE)
  expect_match(err, "ReadImage_fix", fixed = TRUE)
})

test_that("check_sitk_image() gives a generic type error for other wrong types", {
  err <- tryCatch(check_sitk_image(42), error = function(e) conditionMessage(e))
  expect_match(err, "SimpleITK image object")
  expect_match(err, "numeric")
  expect_false(grepl("file path", err))
})

test_that("check_sitk_image() uses the supplied arg_name in the message", {
  err <- tryCatch(check_sitk_image("x.mnc", arg_name = "moving_image"), error = function(e) conditionMessage(e))
  expect_match(err, "`moving_image`", fixed = TRUE)
})

# --- integration: every public function taking `image` rejects a path clearly ---

test_that("public functions taking `image` reject a file path with a clear error", {
  skip_if_not_installed("SimpleITK")

  base <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))
  pkg <- SlicePackage$new(base_slice = base, spacing_k = 1, size_k = 2)
  path <- "brain.nii.gz"

  expect_error(SliceGeometry$from_image_axis(path, "z", 0), "not a file path")
  expect_error(SlicePackage$from_image_axis(path, "z"), "not a file path")
  expect_error(pkg$sample_intensity(path), "not a file path")

  sset <- SlicePackageSet$new(list(only = pkg))
  expect_error(sset$sample_intensity(path), "not a file path")
  expect_error(SlicePackageSet$from_orthogonal_triplet(path, list(x = 0)), "not a file path")

  expect_error(WriteImage_fix(path, "out.nii.gz"), "not a file path")
})
