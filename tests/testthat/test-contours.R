circle_image <- function(n = 20, radius = 5, center = c(10, 10)) {
  filled_sitk_image(c(n, n, n), function(idx) {
    if ((idx[1] - center[1])^2 + (idx[2] - center[2])^2 <= radius^2) 1 else 0
  })
}

two_band_label_image <- function(n = 20) {
  filled_sitk_image(c(n, n, n), function(idx) {
    if (idx[1] < 10) 1 else if (idx[1] < 15) 2 else 0
  })
}

test_that("slice_contours() traces a circle at the expected radius", {
  skip_if_not_installed("SimpleITK")
  img <- circle_image()
  out <- slice_contours(img, axis = "z", coordinate = 10, levels = 0.5)

  expect_true(all(c("package", "k", "level", "obj", "vertex", "x", "y", "z") %in% names(out)))
  expect_true(nrow(out) > 0)
  expect_true(all(out$z == 10))

  dist_from_center <- sqrt((out$x - 10)^2 + (out$y - 10)^2)
  expect_true(all(abs(dist_from_center - 5) < 1.5))
})

test_that("slice_contours() supports multiple levels in one call", {
  skip_if_not_installed("SimpleITK")
  img <- filled_sitk_image(c(20, 20, 20), function(idx) idx[1] / 2)
  out <- slice_contours(img, axis = "z", coordinate = 5, levels = c(2, 5, 8))
  expect_setequal(unique(out$level), c(2, 5, 8))
})

test_that("slice_contours() errors on invalid mask_fill or empty levels", {
  skip_if_not_installed("SimpleITK")
  img <- circle_image()
  expect_error(slice_contours(img, axis = "z", coordinate = 10, levels = 0.5, mask_fill = "nope"))
  expect_error(slice_contours(img, axis = "z", coordinate = 10, levels = numeric(0)))
})

test_that("mask_fill = 'zero' fabricates a boundary at the mask edge; 'nan' does not", {
  skip_if_not_installed("SimpleITK")
  img <- two_band_label_image()
  mask <- filled_sitk_image(c(20, 20, 20), function(idx) if (idx[1] >= 5) 1 else 0)

  out_zero <- slice_contours(img, axis = "z", coordinate = 10, levels = 0.5, mask = mask, mask_fill = "zero")
  out_nan <- slice_contours(img, axis = "z", coordinate = 10, levels = 0.5, mask = mask, mask_fill = "nan")

  # zero-fill: a fabricated edge at the mask boundary (~x=5) PLUS the real edge (~x=15)
  # nan-fill: only the real edge (~x=15)
  expect_true(nrow(out_zero) > nrow(out_nan))
  expect_true(any(out_zero$x < 7))
  expect_false(any(out_nan$x < 7))
})

test_that("slice_label_contours() traces each label's own boundary separately", {
  skip_if_not_installed("SimpleITK")
  img <- two_band_label_image()
  out <- slice_label_contours(img, axis = "z", coordinate = 10)

  expect_setequal(unique(out$label), c(1, 2))
  expect_true(all(c("package", "k", "label", "obj", "vertex", "x", "y", "z") %in% names(out)))
  # label 2 is bounded on both sides (i=10 and i=15) -> more vertices than label 1 (bounded only at i=10)
  expect_gt(sum(out$label == 2), sum(out$label == 1))
})

test_that("slice_label_contours() respects the `labels` subsetting argument", {
  skip_if_not_installed("SimpleITK")
  img <- two_band_label_image()
  out <- slice_label_contours(img, axis = "z", coordinate = 10, labels = 1)
  expect_setequal(unique(out$label), 1)
})

test_that("min_vertices drops small contour paths while keeping larger ones", {
  skip_if_not_installed("SimpleITK")
  # A large circle plus one isolated single-voxel "island" -> a tiny noise contour
  img <- filled_sitk_image(c(30, 30, 30), function(idx) {
    is_big <- (idx[1] - 15)^2 + (idx[2] - 15)^2 <= 64
    is_island <- idx[1] == 2 && idx[2] == 2
    if (is_big || is_island) 1 else 0
  })
  out_unfiltered <- slice_contours(img, axis = "z", coordinate = 15, levels = 0.5)
  n_paths_unfiltered <- nrow(dplyr::distinct(out_unfiltered, package, k, level, obj))
  expect_gt(n_paths_unfiltered, 1)

  out_filtered <- slice_contours(img, axis = "z", coordinate = 15, levels = 0.5, min_vertices = 10)
  n_paths_filtered <- nrow(dplyr::distinct(out_filtered, package, k, level, obj))
  expect_equal(n_paths_filtered, 1)

  path_sizes <- out_filtered %>%
    dplyr::count(package, k, level, obj) %>%
    dplyr::pull(n)
  expect_true(all(path_sizes >= 10))
})

test_that("slice_contours()/slice_label_contours() accept the geometry= escape hatch", {
  skip_if_not_installed("SimpleITK")
  img <- circle_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  out <- slice_contours(img, geometry = geom, levels = 0.5)
  expect_true(nrow(out) > 0)

  label_img <- two_band_label_image()
  geom2 <- SliceGeometry$from_image_axis(label_img, "z", 10)
  outl <- slice_label_contours(label_img, geometry = geom2)
  expect_true(nrow(outl) > 0)
})

test_that("slice_contours()/slice_label_contours() require exactly one of geometry or axis+coordinate", {
  skip_if_not_installed("SimpleITK")
  img <- circle_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  expect_error(slice_contours(img, axis = "z", coordinate = 10, geometry = geom, levels = 0.5))
  expect_error(slice_contours(img, levels = 0.5), "Provide either")
  expect_error(slice_label_contours(img, geometry = geom, axis = "z"))
  expect_error(slice_label_contours(img), "Provide either")
})

test_that("slice_contours() and slice_label_contours() accept file paths for image and mask", {
  skip_if_not_installed("SimpleITK")
  img <- circle_image()
  mask <- filled_sitk_image(c(20, 20, 20), function(idx) 1)
  img_path <- tempfile(fileext = ".nii.gz")
  mask_path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(c(img_path, mask_path)))
  SimpleITK::WriteImage(img, img_path)
  SimpleITK::WriteImage(mask, mask_path)

  out <- slice_contours(img_path, axis = "z", coordinate = 10, levels = 0.5, mask = mask_path)
  expect_true(nrow(out) > 0)
})
