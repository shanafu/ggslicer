# Synthetic MINC-like image: origin/direction chosen to make the x/y negation
# obvious and distinguishable from z (which must never be touched).
make_synthetic_mnc <- function(path, origin = c(-10, -20, -30), value = 42) {
  img <- filled_sitk_image(c(4, 4, 4), function(idx) value, origin = origin)
  SimpleITK::WriteImage(img, path)
}

test_that("orientation_correction() negates x/y of origin and direction, leaves z and pixel data untouched", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(3, 3, 3), function(idx) sum(as.numeric(idx)), origin = c(-10, -20, -30))
  img$SetDirection(c(1, 0, 0, 0, 1, 0, 0, 0, 1))

  corrected <- orientation_correction(img)

  expect_equal(corrected$GetOrigin(), c(10, 20, -30))
  expect_equal(corrected$GetDirection(), c(-1, 0, 0, 0, -1, 0, 0, 0, 1))
  expect_equal(SimpleITK::as.array(corrected), SimpleITK::as.array(img))

  # doesn't mutate the caller's image
  expect_equal(img$GetOrigin(), c(-10, -20, -30))
})

test_that("orientation_correction() is its own inverse", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(3, 3, 3), function(idx) sum(as.numeric(idx)), origin = c(-10, -20, -30))

  twice <- orientation_correction(orientation_correction(img))
  expect_equal(twice$GetOrigin(), img$GetOrigin())
  expect_equal(twice$GetDirection(), img$GetDirection())
  expect_equal(SimpleITK::as.array(twice), SimpleITK::as.array(img))
})

test_that("orientation_correction() also negates a non-identity direction matrix correctly", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(3, 3, 3), function(idx) 0, origin = c(1, 2, 3))
  # a direction matrix with some off-diagonal structure, still orthonormal
  img$SetDirection(c(0, 1, 0, 1, 0, 0, 0, 0, 1))

  corrected <- orientation_correction(img)
  expect_equal(corrected$GetOrigin(), c(-1, -2, 3))
  # flip_diag %*% direction_mat: negate rows 1,2 of the direction matrix (row-major-read)
  expect_equal(corrected$GetDirection(), c(0, -1, 0, -1, 0, 0, 0, 0, 1))
})

test_that("orientation_correction() preserves OriginalFileType metadata when present", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(2, 2, 2), function(idx) 0)
  img$SetMetaData("OriginalFileType", "MINC")
  corrected <- orientation_correction(img)
  expect_equal(corrected$GetMetaData("OriginalFileType"), "MINC")
})

test_that("ReadImage_fix() leaves a non-MINC file's geometry unchanged, tagging it 'Other'", {
  skip_if_not_installed("SimpleITK")
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path))
  make_synthetic_mnc(path) # extension decides format; content doesn't matter here

  raw <- SimpleITK::ReadImage(path)
  fixed <- ReadImage_fix(path)
  expect_equal(fixed$GetOrigin(), raw$GetOrigin())
  expect_equal(fixed$GetDirection(), raw$GetDirection())
  expect_equal(fixed$GetMetaData("OriginalFileType"), "Other")
})

test_that("ReadImage_fix() applies orientation_correction() to a MINC file", {
  skip_if_not_installed("SimpleITK")
  path <- tempfile(fileext = ".mnc")
  on.exit(unlink(path))
  make_synthetic_mnc(path, origin = c(-5, -7, -9))

  raw <- SimpleITK::ReadImage(path)
  fixed <- ReadImage_fix(path)

  expect_equal(fixed$GetOrigin(), c(5, 7, -9))
  expect_equal(fixed$GetMetaData("OriginalFileType"), "MINC")
  expect_equal(SimpleITK::as.array(fixed), SimpleITK::as.array(raw))
})

test_that("WriteImage_fix() round-trips MINC -> MINC back to the exact original raw header", {
  skip_if_not_installed("SimpleITK")
  in_path <- tempfile(fileext = ".mnc")
  out_path <- tempfile(fileext = ".mnc")
  on.exit(unlink(c(in_path, out_path)))
  make_synthetic_mnc(in_path, origin = c(-5, -7, -9))

  raw_original <- SimpleITK::ReadImage(in_path)
  fixed <- ReadImage_fix(in_path)
  WriteImage_fix(fixed, out_path)

  reread <- SimpleITK::ReadImage(out_path)
  expect_equal(reread$GetOrigin(), raw_original$GetOrigin())
  expect_equal(reread$GetDirection(), raw_original$GetDirection())
  expect_equal(SimpleITK::as.array(reread), SimpleITK::as.array(raw_original))
})

test_that("WriteImage_fix() writes MINC -> non-MINC with the corrected header as-is (no restoration)", {
  skip_if_not_installed("SimpleITK")
  in_path <- tempfile(fileext = ".mnc")
  out_path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(c(in_path, out_path)))
  make_synthetic_mnc(in_path, origin = c(-5, -7, -9))

  fixed <- ReadImage_fix(in_path)
  WriteImage_fix(fixed, out_path)

  reread <- SimpleITK::ReadImage(out_path)
  expect_equal(reread$GetOrigin(), fixed$GetOrigin())
  expect_equal(reread$GetDirection(), fixed$GetDirection())
})

test_that("WriteImage_fix() writes an image not originally read as MINC unchanged, regardless of output format", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(3, 3, 3), function(idx) 0, origin = c(1, 2, 3))
  img$SetMetaData("OriginalFileType", "Other")

  out_mnc <- tempfile(fileext = ".mnc")
  on.exit(unlink(out_mnc))
  WriteImage_fix(img, out_mnc)
  reread <- SimpleITK::ReadImage(out_mnc)
  expect_equal(reread$GetOrigin(), img$GetOrigin())
  expect_equal(reread$GetDirection(), img$GetDirection())
})

test_that("WriteImage_fix() errors on a non-SimpleITK-image input", {
  skip_if_not_installed("SimpleITK")
  expect_error(WriteImage_fix("not an image", tempfile(fileext = ".mnc")))
})
