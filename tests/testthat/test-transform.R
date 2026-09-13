test_that("transform_points() applies a translation transform exactly", {
  skip_if_not_installed("SimpleITK")
  df <- tibble::tibble(x = c(0, 1, 2), y = c(0, 0, 0), z = c(0, 0, 0), label = c("a", "b", "c"))
  t <- SimpleITK::TranslationTransform(3L, c(10, 5, -2))
  out <- transform_points(df, t)

  expect_equal(out$x, c(10, 11, 12))
  expect_equal(out$y, c(5, 5, 5))
  expect_equal(out$z, c(-2, -2, -2))
  expect_equal(out$label, df$label) # other columns untouched
})

test_that("transform_points() with invert = TRUE applies the inverse transform", {
  skip_if_not_installed("SimpleITK")
  df <- tibble::tibble(x = 10, y = 5, z = -2)
  t <- SimpleITK::TranslationTransform(3L, c(10, 5, -2))
  out <- transform_points(df, t, invert = TRUE)
  expect_equal(c(out$x, out$y, out$z), c(0, 0, 0))
})

test_that("transform_points() errors clearly when inverting a displacement-field transform", {
  skip_if_not_installed("SimpleITK")
  vec_img <- SimpleITK::Image(as.integer(c(3, 3, 3)), "sitkVectorFloat64", 3L)
  t <- SimpleITK::DisplacementFieldTransform(vec_img)
  df <- tibble::tibble(x = 0, y = 0, z = 0)
  expect_error(transform_points(df, t, invert = TRUE))
})

test_that("transform_points() accepts a transform given as a file path", {
  skip_if_not_installed("SimpleITK")
  t <- SimpleITK::TranslationTransform(3L, c(1, 2, 3))
  path <- tempfile(fileext = ".tfm")
  on.exit(unlink(path))
  SimpleITK::WriteTransform(t, path)

  df <- tibble::tibble(x = 0, y = 0, z = 0)
  out <- transform_points(df, path)
  expect_equal(c(out$x, out$y, out$z), c(1, 2, 3))
})

test_that("transform_points() respects custom x_col/y_col/z_col names", {
  skip_if_not_installed("SimpleITK")
  df <- tibble::tibble(px = 0, py = 0, pz = 0, other = "kept")
  t <- SimpleITK::TranslationTransform(3L, c(1, 2, 3))
  out <- transform_points(df, t, x_col = "px", y_col = "py", z_col = "pz")
  expect_equal(c(out$px, out$py, out$pz), c(1, 2, 3))
  expect_equal(out$other, "kept")
})

# --- read_minc_transform(): synthetic .xfm files, not depending on real testdata ---

write_test_xfm <- function(path, body) {
  writeLines(c("MNI Transform File", "%test", "", body), path)
}

test_that("read_minc_transform() parses a single Linear block", {
  skip_if_not_installed("SimpleITK")
  path <- tempfile(fileext = ".xfm")
  on.exit(unlink(path))
  write_test_xfm(path, c(
    "Transform_Type = Linear;",
    "Linear_Transform =",
    " 1 0 0 5",
    " 0 1 0 6",
    " 0 0 1 7;"
  ))

  t <- read_minc_transform(path)
  expect_true(inherits(t, "_p_itk__simple__Transform"))
  out <- t$TransformPoint(c(0, 0, 0))
  expect_equal(out, c(5, 6, 7))
})

test_that("read_minc_transform() parses a single Grid_Transform block", {
  skip_if_not_installed("SimpleITK")
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  grid_path <- file.path(tmpdir, "grid.mnc")
  vec_img <- SimpleITK::Image(as.integer(c(3, 3, 3)), "sitkVectorFloat64", 3L)
  for (i in 0:2) for (j in 0:2) for (k in 0:2) {
    vec_img$SetPixel(c(i, j, k), c(1, 2, 3))
  }
  SimpleITK::WriteImage(vec_img, grid_path)

  xfm_path <- file.path(tmpdir, "test.xfm")
  write_test_xfm(xfm_path, c(
    "Transform_Type = Grid_Transform;",
    "Displacement_Volume = grid.mnc;"
  ))

  t <- read_minc_transform(xfm_path)
  out <- t$TransformPoint(c(0, 0, 0))
  expect_equal(out, c(1, 2, 3))
})

test_that("read_minc_transform() concatenates multiple blocks in the correct order", {
  skip_if_not_installed("SimpleITK")
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  grid_path <- file.path(tmpdir, "grid.mnc")
  vec_img <- SimpleITK::Image(as.integer(c(3, 3, 3)), "sitkVectorFloat64", 3L)
  for (i in 0:2) for (j in 0:2) for (k in 0:2) {
    vec_img$SetPixel(c(i, j, k), c(100, 0, 0))
  }
  SimpleITK::WriteImage(vec_img, grid_path)

  xfm_path <- file.path(tmpdir, "test.xfm")
  write_test_xfm(xfm_path, c(
    "Transform_Type = Linear;",
    "Linear_Transform =",
    " 2 0 0 0",
    " 0 2 0 0",
    " 0 0 2 0;",
    "Transform_Type = Grid_Transform;",
    "Displacement_Volume = grid.mnc;"
  ))

  t <- read_minc_transform(xfm_path)
  expect_true(inherits(t, "_p_itk__simple__CompositeTransform"))

  # file order [Linear, Grid]: apply linear (scale by 2) first, then grid (+100 in x)
  out <- t$TransformPoint(c(1, 0, 0))
  expect_equal(out, c(102, 0, 0))
})

test_that("read_minc_transform() resolves Displacement_Volume relative to the .xfm's own directory", {
  skip_if_not_installed("SimpleITK")
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  vec_img <- SimpleITK::Image(as.integer(c(2, 2, 2)), "sitkVectorFloat64", 3L)
  SimpleITK::WriteImage(vec_img, file.path(tmpdir, "somegrid.mnc"))

  xfm_path <- file.path(tmpdir, "test.xfm")
  write_test_xfm(xfm_path, c(
    "Transform_Type = Grid_Transform;",
    "Displacement_Volume = somegrid.mnc;"
  ))

  expect_no_error(read_minc_transform(xfm_path))
})

test_that("read_minc_transform() errors clearly on malformed or unsupported input", {
  skip_if_not_installed("SimpleITK")
  not_an_xfm <- tempfile(fileext = ".xfm")
  on.exit(unlink(not_an_xfm), add = TRUE)
  writeLines("not a transform file", not_an_xfm)
  expect_error(read_minc_transform(not_an_xfm))

  bad_linear <- tempfile(fileext = ".xfm")
  on.exit(unlink(bad_linear), add = TRUE)
  write_test_xfm(bad_linear, c("Transform_Type = Linear;", "Linear_Transform =", " 1 2 3;"))
  expect_error(read_minc_transform(bad_linear))

  unsupported <- tempfile(fileext = ".xfm")
  on.exit(unlink(unsupported), add = TRUE)
  write_test_xfm(unsupported, "Transform_Type = Thin_Plate_Spline_Transform;")
  expect_error(read_minc_transform(unsupported), "Unsupported")
})

test_that("transform_points() routes .xfm paths through read_minc_transform()", {
  skip_if_not_installed("SimpleITK")
  path <- tempfile(fileext = ".xfm")
  on.exit(unlink(path))
  write_test_xfm(path, c(
    "Transform_Type = Linear;",
    "Linear_Transform =",
    " 1 0 0 1",
    " 0 1 0 2",
    " 0 0 1 3;"
  ))

  df <- tibble::tibble(x = 0, y = 0, z = 0)
  out <- transform_points(df, path)
  expect_equal(c(out$x, out$y, out$z), c(1, 2, 3))
})
