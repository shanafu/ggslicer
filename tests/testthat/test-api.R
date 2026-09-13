test_that("discrete_data_names() includes the agreed privileged names", {
  nms <- discrete_data_names()
  expect_true(all(c("mask", "label", "labels", "segmentation", "segmentations", "atlas") %in% nms))
})

test_that(".resolve_interpolator() matches whole tokens, not substrings", {
  expect_equal(.resolve_interpolator("brain_mask", "sitkLinear", discrete_data_names()), "sitkNearestNeighbor")
  expect_equal(.resolve_interpolator("aseg.mgz", "sitkLinear", discrete_data_names()), "sitkNearestNeighbor")
  expect_equal(.resolve_interpolator("landmasking_score", "sitkLinear", discrete_data_names()), "sitkLinear")
  expect_equal(.resolve_interpolator("tstat1", "sitkLinear", discrete_data_names()), "sitkLinear")
})

test_that("sample_images() defaults a single non-list image to an `intensity` column", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) idx[1] + 10 * idx[2] + 100 * idx[3])
  geom <- SliceGeometry$from_image_axis(img, "z", 1)
  out <- sample_images(geom, img)
  expect_true("intensity" %in% names(out))
})

test_that("sample_images() names columns after the given list names", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) idx[1] + 10 * idx[2] + 100 * idx[3])
  mask <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric((idx[1] + idx[2]) %% 2 == 0))
  geom <- SliceGeometry$from_image_axis(img, "z", 1)
  out <- sample_images(geom, list(value = img, mask = mask))
  expect_true(all(c("value", "mask") %in% names(out)))
  expect_false("intensity" %in% names(out))
})

test_that("sample_images() accepts a file path and reads it internally", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) idx[1] + 10 * idx[2] + 100 * idx[3])
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path))
  SimpleITK::WriteImage(img, path)
  geom <- SliceGeometry$from_image_axis(img, "z", 1)
  out <- sample_images(geom, path)
  expect_true("intensity" %in% names(out))
  expect_equal(nrow(out), 16)
})

test_that("interpolator_overrides wins over the name-based rule", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) idx[1] + 10 * idx[2] + 100 * idx[3])
  geom <- SliceGeometry$from_image_axis(img, "z", 1)
  expect_no_error(
    sample_images(geom, list(value = img), interpolator_overrides = list(value = "sitkNearestNeighbor"))
  )
})

test_that("sample_images() broadcasts a 3D image's values across extra_index combinations from a 4D image", {
  skip_if_not_installed("SimpleITK")
  main4d <- filled_sitk_image(c(4, 4, 4, 2), function(idx) idx[1] + 10 * idx[2] + 100 * idx[3] + 1000 * idx[4])
  mask3d <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric((idx[1] + idx[2]) %% 2 == 0))
  geom <- SliceGeometry$from_image_axis(mask3d, "z", 1)
  out <- sample_images(
    geom, list(value = main4d, mask = mask3d),
    extra_index = list(t = 0:1)
  )
  expect_true(all(c("t", "value", "mask") %in% names(out)))
  expect_equal(nrow(out), 16 * 2)
  by_point <- split(out$mask, paste(out$i, out$j, out$k))
  expect_true(all(vapply(by_point, function(v) length(unique(v)) == 1, logical(1))))
})

test_that("sample_images() errors clearly when an image's dimensionality doesn't match extra_index", {
  skip_if_not_installed("SimpleITK")
  ref3d <- filled_sitk_image(c(2, 2, 2), function(idx) as.numeric(sum(idx)))
  main5d <- filled_sitk_image(c(2, 2, 2, 2, 2), function(idx) as.numeric(sum(idx)))
  geom <- SliceGeometry$from_image_axis(ref3d, "z", 0)
  expect_error(
    sample_images(geom, list(value = main5d), extra_index = list(t = 0:1)),
    "extra_index"
  )
})

test_that("build_slice_geometry() returns one package for evenly spaced coordinates", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric(sum(idx)))
  set <- build_slice_geometry(img, "z", c(0, 1, 2))
  expect_s3_class(set, "SlicePackageSet")
  expect_length(set$get_package_names(), 1)
})

test_that("build_slice_geometry() returns one package per slice for unevenly spaced coordinates", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric(sum(idx)))
  set <- build_slice_geometry(img, "z", c(0, 1, 3))
  expect_length(set$get_package_names(), 3)
})

test_that("build_slice_geometry() with a single coordinate returns a single-slice package", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric(sum(idx)))
  set <- build_slice_geometry(img, "z", 1)
  expect_length(set$get_package_names(), 1)
})

test_that("slice_image() requires exactly one of geometry or axis+coordinate", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric(sum(idx)))
  geom <- SliceGeometry$from_image_axis(img, "z", 1)
  expect_error(slice_image(img, axis = "z", coordinate = 1, geometry = geom))
  expect_error(slice_image(img), "Provide either")
  expect_error(slice_image(img, axis = "z"), "Provide either")
})

test_that("slice_image() samples the main image plus extra images together", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) idx[1] + 10 * idx[2] + 100 * idx[3])
  mask <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric((idx[1] + idx[2]) %% 2 == 0))
  out <- slice_image(img, axis = "z", coordinate = 1, extra_images = list(mask = mask))
  expect_true(all(c("value", "mask", "package", "x", "y", "z") %in% names(out)))
})

test_that("slice_image() accepts the geometry= escape hatch", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(4, 4, 4), function(idx) as.numeric(sum(idx)))
  geom <- SliceGeometry$from_image_axis(img, "z", 1)
  out <- slice_image(img, geometry = geom)
  expect_true("value" %in% names(out))
})
