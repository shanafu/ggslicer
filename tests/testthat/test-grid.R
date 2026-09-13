simple_image <- function(n = 20) {
  filled_sitk_image(c(n, n, n), function(idx) 0)
}

test_that("slice_grid() returns the expected columns and both parts", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  out <- slice_grid(geometry = geom)

  expect_true(all(c("package", "k", "part", "grid_axis", "line_id", "vertex", "x", "y", "z") %in% names(out)))
  expect_setequal(unique(out$part), c("box", "grid"))
  expect_setequal(unique(out$grid_axis), c("i", "j"))
})

test_that("the box sits exactly at the slice's true extent when padding = 0", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  extent <- geom$get_extent()

  out <- slice_grid(geometry = geom, padding = 0)
  box <- out[out$part == "box", ]

  i_lines <- box[box$grid_axis == "i", ]
  i_coords <- sort(unique(round(sapply(split(i_lines$x, i_lines$line_id), function(v) v[1]), 6)))
  # For an axial slice (in_plane axes = x,y with direction_i = x, direction_j = y here),
  # the i-line grid coordinate is reflected directly in the constant world x value.
  expect_equal(length(i_coords), 2)
  expect_equal(diff(range(box$x)), extent[1], tolerance = 1e-6)
  expect_equal(diff(range(box$y)), extent[2], tolerance = 1e-6)
})

test_that("padding extends the box beyond the slice's true extent", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  extent <- geom$get_extent()

  out <- slice_grid(geometry = geom, padding = 2)
  box <- out[out$part == "box", ]
  expect_equal(diff(range(box$x)), extent[1] + 4, tolerance = 1e-6)
  expect_equal(diff(range(box$y)), extent[2] + 4, tolerance = 1e-6)
})

test_that("default spacing is scale-adaptive: min(extent)/10", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image(40)
  geom <- SliceGeometry$from_image_axis(img, "z", 20)
  extent <- geom$get_extent()
  expected_spacing <- min(extent) / 10

  out <- slice_grid(geometry = geom)
  grid_only <- out[out$part == "grid" & out$grid_axis == "i", ]
  grid_coords <- sort(unique(round(sapply(split(grid_only$y, grid_only$line_id), function(v) v[1]), 6)))
  observed_spacing <- diff(grid_coords)
  expect_true(all(abs(observed_spacing - expected_spacing) < 1e-6))
})

test_that("explicit spacing and point_spacing control line count and point density", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image(40)
  geom <- SliceGeometry$from_image_axis(img, "z", 20)
  extent <- geom$get_extent()

  out <- slice_grid(geometry = geom, spacing = 5, point_spacing = 1, padding = 0)
  grid_i <- out[out$part == "grid" & out$grid_axis == "i", ]
  n_lines <- length(unique(grid_i$line_id))
  expect_equal(n_lines, length(seq(0, extent[1], by = 5)))

  one_line <- grid_i[grid_i$line_id == 1, ]
  expect_equal(nrow(one_line), round(extent[2] / 1) + 1)
})

test_that("each slice in a multi-slice SlicePackage gets its own grid", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  pkg <- SlicePackage$from_image_axis(img, "z")
  pkg$set_size_k(3)
  out <- slice_grid(geometry = pkg)
  expect_setequal(unique(out$k), 0:2)
  # each k should have the same number of rows (identical geometry per slice)
  counts <- table(out$k)
  expect_true(length(unique(as.numeric(counts))) == 1)
})

test_that("each package in a SlicePackageSet gets its own independent grid", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  set <- SlicePackageSet$from_orthogonal_triplet(img, list(x = 0, y = 0, z = 0))
  out <- slice_grid(geometry = set)
  expect_setequal(unique(out$package), c("sagittal", "coronal", "axial"))
})

test_that("slice_grid() requires exactly one of geometry or axis+coordinate", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  expect_error(slice_grid(img, axis = "z", coordinate = 10, geometry = geom))
  expect_error(slice_grid(), "Provide either")
  expect_error(slice_grid(axis = "z"), "Provide either")
})

test_that("slice_grid() requires image when geometry is not supplied", {
  skip_if_not_installed("SimpleITK")
  expect_error(slice_grid(axis = "z", coordinate = 10), "image.*required")
})

test_that("slice_grid() validates spacing/point_spacing/padding", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  expect_error(slice_grid(geometry = geom, spacing = -1))
  expect_error(slice_grid(geometry = geom, point_spacing = 0))
  expect_error(slice_grid(geometry = geom, padding = -1000))
})

test_that("slice_grid_layers() returns two styled geom_path layers that build into a plot", {
  skip_if_not_installed("SimpleITK")
  img <- simple_image()
  geom <- SliceGeometry$from_image_axis(img, "z", 10)
  grid_df <- slice_grid(geometry = geom)

  layers <- slice_grid_layers(grid_df, box_color = "blue", grid_color = "green")
  expect_true(all(c("box", "grid") %in% names(layers)))
  expect_s3_class(layers$box, "ggproto")
  expect_s3_class(layers$grid, "ggproto")
  expect_equal(layers$box$aes_params$colour, "blue")
  expect_equal(layers$grid$aes_params$colour, "green")

  plt <- ggplot2::ggplot(grid_df, ggplot2::aes(x = x, y = y)) + layers$box + layers$grid
  expect_s3_class(plt, "ggplot")
})
