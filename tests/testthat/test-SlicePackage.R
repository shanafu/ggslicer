base_slice <- function(...) {
  defaults <- list(origin = c(0, 0, 0), direction_i = c(1, 0, 0), direction_j = c(0, 1, 0), spacing = c(1, 1), size = c(3, 3))
  args <- utils::modifyList(defaults, list(...))
  do.call(SliceGeometry$new, args)
}

test_that("constructor validates its inputs", {
  expect_error(SlicePackage$new(base_slice = "not a slice", spacing_k = 1, size_k = 2), "SliceGeometry")
  expect_error(SlicePackage$new(base_slice = base_slice(), spacing_k = 0, size_k = 2), "> 0")
  expect_error(SlicePackage$new(base_slice = base_slice(), spacing_k = -1, size_k = 2), "> 0")
  expect_error(SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 0), ">= 1")
  expect_error(SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 1.5), ">= 1")
})

test_that("size/spacing/normal combine the base slice with the k axis", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 2, size_k = 4)

  expect_equal(pkg$get_size(), c(3L, 3L, 4L))
  expect_equal(pkg$get_spacing(), c(1, 1, 2))
  expect_equal(unname(pkg$get_normal()), c(0, 0, 1))
})

test_that("get_slice() returns the correctly-translated k-th slice", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 2, size_k = 4)

  expect_equal(pkg$get_slice(0)$get_origin(), c(0, 0, 0))
  expect_equal(pkg$get_slice(2)$get_origin(), c(0, 0, 4))
  expect_equal(pkg$get_slice(2)$get_size(), pkg$get_base_slice()$get_size())

  expect_error(pkg$get_slice(-1), "0:\\(size_k")
  expect_error(pkg$get_slice(4), "0:\\(size_k")
  expect_error(pkg$get_slice(1.5), "0:\\(size_k")
})

test_that("get_sample_points() produces the full 3D stack", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 2, size_k = 4)
  pts <- pkg$get_sample_points()

  expect_equal(nrow(pts), 3 * 3 * 4)
  expect_named(pts, c("i", "j", "k", "x", "y", "z"))

  pt <- pts[pts$i == 0 & pts$j == 0 & pts$k == 2, c("x", "y", "z")]
  expect_equal(unname(unlist(pt)), c(0, 0, 4))
})

test_that("get_sample_points() is memoized and set_spacing_k()/set_size_k() invalidate it", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 2)
  pts1 <- pkg$get_sample_points()
  expect_identical(pts1, pkg$get_sample_points())

  pkg$set_spacing_k(3)
  expect_equal(max(pkg$get_sample_points()$z), 3)

  pkg$set_size_k(5)
  expect_equal(nrow(pkg$get_sample_points()), 3 * 3 * 5)
})

test_that("as_sitk_reference_image() spans the whole stack in one geometry", {
  skip_if_not_installed("SimpleITK")

  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 2, size_k = 4)
  ref <- pkg$as_sitk_reference_image()

  expect_equal(ref$GetSize(), c(3L, 3L, 4L))
  expect_equal(ref$GetSpacing(), c(1, 1, 2))
  expect_equal(ref$GetDirection(), c(1, 0, 0, 0, 1, 0, 0, 0, 1))
})

test_that("sample_intensity() resamples a 3D image in one call and matches source pixels", {
  skip_if_not_installed("SimpleITK")

  pkg <- SlicePackage$new(base_slice = base_slice(size = c(5, 5)), spacing_k = 1, size_k = 10)
  src <- filled_sitk_image(c(5, 5, 10), function(idx) 100 * idx[3] + 10 * idx[1] + idx[2])

  res <- pkg$sample_intensity(src)
  expect_named(res, c("i", "j", "k", "x", "y", "z", "intensity"))

  row <- res[res$i == 2 & res$j == 1 & res$k == 1, ]
  expect_equal(row$intensity, 100 * 1 + 10 * 2 + 1)
})

test_that("sample_intensity() marks out-of-bounds samples as NA, not NaN", {
  skip_if_not_installed("SimpleITK")

  small_src <- filled_sitk_image(c(5, 5, 10), function(idx) 1)
  big_pkg <- SlicePackage$new(base_slice = base_slice(size = c(10, 10)), spacing_k = 1, size_k = 1)

  res <- big_pkg$sample_intensity(small_src)
  expect_true(any(is.na(res$intensity)))
  expect_false(any(is.nan(res$intensity)))
})

test_that("sample_intensity() rejects non-3D images", {
  skip_if_not_installed("SimpleITK")

  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 2)
  img4d <- filled_sitk_image(c(2, 2, 2, 2), function(idx) 0)

  expect_error(pkg$sample_intensity(img4d), "3D image")
})

test_that("get_base_slice() returns a clone, not the live internal slice", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 2)
  pkg$get_sample_points() # populate the cache

  clone <- pkg$get_base_slice()
  clone$set_origin(c(99, 99, 99))

  expect_equal(pkg$get_base_slice()$get_origin(), c(0, 0, 0))
  expect_equal(pkg$get_sample_points()$x[1], 0) # cache untouched by the external mutation
})

test_that("set_base_slice() replaces the base slice and invalidates the cache", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 2)
  pkg$get_sample_points()

  new_base <- SliceGeometry$new(c(1, 1, 1), c(1, 0, 0), c(0, 1, 0), c(2, 2), c(4, 4))
  pkg$set_base_slice(new_base)

  expect_equal(pkg$get_base_slice()$get_origin(), c(1, 1, 1))
  expect_equal(pkg$get_size(), c(4L, 4L, 2L))
  expect_error(pkg$set_base_slice("not a slice"), "SliceGeometry")
})

test_that("set_spacing()/set_size() update all 3 axes at once", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 2)

  pkg$set_spacing(c(2, 3, 4))
  expect_equal(pkg$get_spacing(), c(2, 3, 4))

  pkg$set_size(c(5, 6, 7))
  expect_equal(pkg$get_size(), c(5L, 6L, 7L))

  expect_error(pkg$set_spacing(c(1, 1)), "length-3")
  expect_error(pkg$set_size(c(1, 1)), "length-3")
})

test_that("SlicePackage$from_extent_k() derives spacing_k from a total thickness", {
  pkg <- SlicePackage$from_extent_k(base_slice = base_slice(), extent_k = 10, size_k = 6)
  expect_equal(pkg$get_spacing_k(), 2)
  expect_equal(pkg$get_size_k(), 6L)

  expect_error(SlicePackage$from_extent_k(base_slice(), extent_k = 10, size_k = 1), ">= 2")
  expect_error(SlicePackage$from_extent_k(base_slice(), extent_k = -1, size_k = 6), "> 0")
})

test_that("SlicePackage$from_center_k() treats the given slice as the stack's middle", {
  center <- base_slice()
  pkg <- SlicePackage$from_center_k(center_slice = center, spacing_k = 2, size_k = 5)

  middle <- pkg$get_slice(2) # 0-indexed: k=2 is the middle of 5 slices
  expect_equal(middle$get_origin(), center$get_origin())
  expect_equal(pkg$get_slice(0)$get_origin(), center$get_origin() - 4 * unname(center$get_normal()))
})

test_that("SlicePackage$print() reports size, spacing, normal, and base-slice geometry", {
  pkg <- SlicePackage$new(base_slice = base_slice(), spacing_k = 1, size_k = 2)
  expect_output(print(pkg), "<SlicePackage>")
  expect_output(print(pkg), "normal:")
  expect_output(print(pkg), "base origin:")
})

test_that("SlicePackage$from_slices() adopts a list of pre-built parallel slices", {
  base <- base_slice()
  slices <- lapply(0:3, function(k) base$translate(k * 2 * base$get_normal()))

  pkg <- SlicePackage$from_slices(slices)
  expect_equal(pkg$get_spacing_k(), 2)
  expect_equal(pkg$get_size_k(), 4L)
  expect_equal(pkg$get_base_slice()$get_origin(), base$get_origin())

  expect_error(SlicePackage$from_slices(list(base)), "at least 2")
  expect_error(SlicePackage$from_slices(list(base, "not a slice")), "SliceGeometry")

  different_size <- base$clone()
  different_size$set_size(c(5, 5))
  expect_error(SlicePackage$from_slices(list(base, different_size)), "same size")

  uneven <- slices
  uneven[[3]] <- uneven[[3]]$translate(c(0, 0, 0.5))
  expect_error(SlicePackage$from_slices(uneven), "evenly spaced")

  off_axis <- slices
  off_axis[[2]] <- off_axis[[2]]$translate(c(0.1, 0, 0))
  expect_error(SlicePackage$from_slices(off_axis), "aligned along the shared normal")
})

test_that("SlicePackage$from_image_axis() spans the whole image at its own resolution", {
  skip_if_not_installed("SimpleITK")

  img <- filled_sitk_image(c(6, 7, 8), function(idx) 0, origin = c(-3, -3.5, -4), spacing = c(1, 1, 1))
  pkg <- SlicePackage$from_image_axis(img, "axial")

  expect_equal(pkg$get_size(), c(6L, 7L, 8L))
  expect_equal(pkg$get_spacing(), c(1, 1, 1))
  expect_equal(pkg$get_base_slice()$get_origin(), img$GetOrigin())

  pkg_num <- SlicePackage$from_image_axis(img, 3)
  expect_equal(pkg_num$get_size(), pkg$get_size())

  img4d <- filled_sitk_image(c(2, 2, 2, 2), function(idx) 0)
  expect_error(SlicePackage$from_image_axis(img4d, "z"), "3D image")
})
