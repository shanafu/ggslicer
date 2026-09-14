make_vector_image <- function(size, value_fn, origin = NULL, spacing = NULL) {
  vec_img <- SimpleITK::Image(as.integer(size), "sitkVectorFloat64", 3L)
  if (!is.null(origin)) vec_img$SetOrigin(origin)
  if (!is.null(spacing)) vec_img$SetSpacing(spacing)
  for (i in 0:(size[1] - 1)) for (j in 0:(size[2] - 1)) for (k in 0:(size[3] - 1)) {
    vec_img$SetPixel(c(i, j, k), value_fn(c(i, j, k)))
  }
  vec_img
}

write_test_xfm <- function(path, body) {
  writeLines(c("MNI Transform File", "%test", "", body), path)
}

test_that("a constant translation warp gives the expected in-plane/normal decomposition on an axial slice", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(20, 20, 5), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(3, 4, 5))

  arrows <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 5)

  expect_true(all(c(
    "package", "k", "i", "j", "x", "y", "z", "xend", "yend", "zend",
    "in_plane_displacement", "normal_displacement", "total_displacement"
  ) %in% names(arrows)))

  # axial slice: direction_i = x, direction_j = y, normal = z (up to sign)
  expect_equal(arrows$in_plane_displacement, rep(5, nrow(arrows))) # sqrt(3^2+4^2) = 5
  expect_equal(abs(arrows$normal_displacement), rep(5, nrow(arrows)))
  expect_equal(arrows$total_displacement, sqrt(arrows$in_plane_displacement^2 + arrows$normal_displacement^2))

  expect_equal(arrows$xend - arrows$x, rep(3, nrow(arrows)))
  expect_equal(arrows$yend - arrows$y, rep(4, nrow(arrows)))
  expect_equal(arrows$zend, arrows$z) # arrow never leaves the slice plane
})

test_that("arrow_length rescales the drawn arrow but not the reported true displacement", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(20, 20, 5), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(3, 4, 0))

  true_arrows <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 5)
  scaled_arrows <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 5, arrow_length = 2)

  expect_equal(true_arrows$in_plane_displacement, scaled_arrows$in_plane_displacement)
  expect_equal(true_arrows$total_displacement, scaled_arrows$total_displacement)

  drawn_len <- sqrt(
    (scaled_arrows$xend - scaled_arrows$x)^2 +
      (scaled_arrows$yend - scaled_arrows$y)^2 +
      (scaled_arrows$zend - scaled_arrows$z)^2
  )
  expect_equal(drawn_len, rep(2, nrow(scaled_arrows)))

  # direction preserved: unit vector of the drawn arrow matches the true displacement's
  true_dir <- cbind(true_arrows$xend - true_arrows$x, true_arrows$yend - true_arrows$y)
  true_dir <- true_dir / sqrt(rowSums(true_dir^2))
  scaled_dir <- cbind(scaled_arrows$xend - scaled_arrows$x, scaled_arrows$yend - scaled_arrows$y)
  scaled_dir <- scaled_dir / sqrt(rowSums(scaled_dir^2))
  expect_equal(scaled_dir, true_dir)
})

test_that("a purely out-of-plane displacement leaves a zero-length arrow even with arrow_length set", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(20, 20, 5), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(0, 0, 7)) # pure z (normal) displacement on an axial slice

  arrows <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 5, arrow_length = 3)

  expect_equal(arrows$in_plane_displacement, rep(0, nrow(arrows)))
  expect_equal(arrows$normal_displacement, rep(7, nrow(arrows)))
  expect_equal(arrows$xend, arrows$x)
  expect_equal(arrows$yend, arrows$y)
  expect_equal(arrows$zend, arrows$z)
})

test_that("a raw vector image and an already-built DisplacementFieldTransform give identical results", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(10, 10, 3), function(idx) 0)
  vec_img_for_image_form <- make_vector_image(c(10, 10, 3), function(idx) c(idx[1] * 0.1, idx[2] * 0.2, 0.5))
  arrows_image_form <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = vec_img_for_image_form, spacing = 3)

  # DisplacementFieldTransform(image) moves/invalidates the source image object
  # (confirmed in this codebase; see warpfield.R and CLAUDE.md) -- build a fresh
  # copy for the Transform-object input form.
  vec_img_for_transform_form <- make_vector_image(c(10, 10, 3), function(idx) c(idx[1] * 0.1, idx[2] * 0.2, 0.5))
  transform_obj <- SimpleITK::DisplacementFieldTransform(vec_img_for_transform_form)
  arrows_transform_form <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = transform_obj, spacing = 3)

  expect_equal(arrows_image_form$total_displacement, arrows_transform_form$total_displacement)
  expect_equal(arrows_image_form$xend, arrows_transform_form$xend)
})

test_that("warp given as a .xfm path (Linear block) matches the equivalent Transform object", {
  skip_if_not_installed("SimpleITK")
  path <- tempfile(fileext = ".xfm")
  on.exit(unlink(path))
  write_test_xfm(path, c(
    "Transform_Type = Linear;",
    "Linear_Transform =",
    " 1 0 0 2",
    " 0 1 0 -1",
    " 0 0 1 0;"
  ))

  image <- filled_sitk_image(c(10, 10, 3), function(idx) 0)
  arrows_path <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = path, spacing = 3)

  t <- SimpleITK::TranslationTransform(3L, c(2, -1, 0))
  arrows_transform <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 3)

  expect_equal(arrows_path$xend, arrows_transform$xend)
  expect_equal(arrows_path$yend, arrows_transform$yend)
})

test_that("invert = TRUE applies the inverse warp", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(10, 10, 3), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(3, 4, 0))

  forward <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 5)
  inverted <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 5, invert = TRUE)

  expect_equal(inverted$xend - inverted$x, rep(-3, nrow(inverted)))
  expect_equal(inverted$yend - inverted$y, rep(-4, nrow(inverted)))
  expect_equal(forward$in_plane_displacement, inverted$in_plane_displacement)
})

test_that("spacing controls lattice density and defaults to slice_grid()'s own min(extent)/10 convention", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(40, 40, 3), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(1, 1, 0))

  coarse <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 10)
  fine <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 2)
  expect_lt(nrow(coarse), nrow(fine))

  geom <- SliceGeometry$from_image_axis(image, "axial", 0)
  default_arrows <- slice_warp_arrows(geometry = geom, warp = t)
  extent <- geom$get_extent()
  expected_spacing <- min(extent) / 10
  expected_size <- round(extent / expected_spacing) + 1
  expect_equal(nrow(default_arrows), prod(expected_size))
})

test_that("slice_warp_arrows() supports a multi-slice SlicePackage and the geometry= escape hatch", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(10, 10, 5), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(1, 1, 1))
  pkg <- SlicePackage$from_image_axis(image, "axial")

  arrows <- slice_warp_arrows(geometry = pkg, warp = t, spacing = 3)
  expect_setequal(unique(arrows$k), 0:4)
  expect_true(all(arrows$package == "package_1"))
})

test_that("slice_warp_arrows() validates its arguments", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(10, 10, 3), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(1, 1, 1))

  expect_error(slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = -1))
  expect_error(slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, arrow_length = 0))
  expect_error(slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = "not_a_real_object"))
  expect_error(slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, geometry = SliceGeometry$from_image_axis(image, "axial", 0)))
})

test_that("the in-plane/normal decomposition satisfies the orthonormal invariant on an oblique slice", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(20, 20, 20), function(idx) 0)
  vec_img <- make_vector_image(c(20, 20, 20), function(idx) c(idx[1] * 0.05, -idx[2] * 0.03, idx[3] * 0.02))
  geom <- SliceGeometry$from_normal(
    origin = c(2, 2, 2), normal = c(1, 1, 1), spacing = c(1, 1), size = c(5, 5)
  )

  arrows <- slice_warp_arrows(geometry = geom, warp = vec_img, spacing = 1)
  inv_err <- max(abs(
    arrows$total_displacement^2 -
      (arrows$in_plane_displacement^2 + arrows$normal_displacement^2)
  ))
  expect_lt(inv_err, 1e-8)
})

test_that("slice_warp_arrows_layer() builds a usable ggplot2 layer", {
  skip_if_not_installed("SimpleITK")
  image <- filled_sitk_image(c(10, 10, 3), function(idx) 0)
  t <- SimpleITK::TranslationTransform(3L, c(1, 1, 0))
  arrows <- slice_warp_arrows(image, axis = "axial", coordinate = 0, warp = t, spacing = 3)

  layer <- slice_warp_arrows_layer(arrows)
  expect_s3_class(layer, "ggproto")

  p <- ggplot2::ggplot(arrows, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) + layer
  expect_s3_class(p, "ggplot")
})
