test_that("primary constructor validates its inputs", {
  ok <- function(...) SliceGeometry$new(...)

  expect_error(
    ok(origin = c(0, 0), direction_i = c(1, 0, 0), direction_j = c(0, 1, 0), spacing = c(1, 1), size = c(2, 2)),
    "length-3"
  )
  expect_error(
    ok(origin = c(0, 0, 0), direction_i = c(2, 0, 0), direction_j = c(0, 1, 0), spacing = c(1, 1), size = c(2, 2)),
    "unit vector"
  )
  expect_error(
    ok(origin = c(0, 0, 0), direction_i = c(1, 0, 0), direction_j = c(1, 0, 0), spacing = c(1, 1), size = c(2, 2)),
    "orthogonal"
  )
  expect_error(
    ok(origin = c(0, 0, 0), direction_i = c(1, 0, 0), direction_j = c(0, 1, 0), spacing = c(0, 1), size = c(2, 2)),
    "spacing"
  )
  expect_error(
    ok(origin = c(0, 0, 0), direction_i = c(1, 0, 0), direction_j = c(0, 1, 0), spacing = c(1, 1), size = c(0, 2)),
    "size"
  )
  expect_error(
    ok(origin = c(0, 0, 0), direction_i = c(1, 0, 0), direction_j = c(0, 1, 0), spacing = c(1, 1), size = c(2.5, 2)),
    "size"
  )
})

test_that("getters return the geometry that was constructed", {
  s <- SliceGeometry$new(
    origin = c(1, 2, 3),
    direction_i = c(1, 0, 0),
    direction_j = c(0, 1, 0),
    spacing = c(1, 2),
    size = c(3, 4)
  )

  expect_equal(s$get_origin(), c(1, 2, 3))
  expect_equal(unname(s$get_direction_i()), c(1, 0, 0))
  expect_equal(unname(s$get_direction_j()), c(0, 1, 0))
  expect_equal(s$get_spacing(), c(1, 2))
  expect_equal(s$get_size(), c(3L, 4L))
  expect_equal(unname(s$get_normal()), c(0, 0, 1))
  expect_equal(s$get_extent(), c(2, 6))
  expect_equal(unname(s$get_center()), c(1 + 1, 2 + 3, 3))
})

test_that("get_bounds() returns the four corners", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 2), c(3, 4))
  b <- s$get_bounds()

  expect_equal(nrow(b), 4)
  expect_setequal(names(b), c("corner", "x", "y", "z"))
  expect_equal(
    as.matrix(b[b$corner == "i1_j1", c("x", "y", "z")]),
    matrix(c(2, 6, 0), nrow = 1, dimnames = list(NULL, c("x", "y", "z")))
  )
})

test_that("get_sample_points() is correct and 0-indexed", {
  s <- SliceGeometry$new(
    origin = c(0, 0, 0),
    direction_i = c(1, 0, 0),
    direction_j = c(0, 1, 0),
    spacing = c(1, 2),
    size = c(3, 4)
  )
  pts <- s$get_sample_points()

  expect_equal(nrow(pts), 12)
  expect_named(pts, c("i", "j", "x", "y", "z"))
  expect_equal(min(pts$i), 0)
  expect_equal(max(pts$i), 2)
  expect_equal(min(pts$j), 0)
  expect_equal(max(pts$j), 3)

  corner <- pts[pts$i == 2 & pts$j == 3, c("x", "y", "z")]
  expect_equal(unname(unlist(corner)), c(2, 6, 0))
})

test_that("get_sample_points() handles an oblique, non-axis-aligned plane", {
  # direction_i / direction_j chosen as an orthonormal pair not aligned to any
  # cartesian axis
  di <- c(1, 1, 0) / sqrt(2)
  dj <- c(-1, 1, 0) / sqrt(2)
  s <- SliceGeometry$new(origin = c(5, 5, 5), direction_i = di, direction_j = dj, spacing = c(1, 1), size = c(2, 2))
  pts <- s$get_sample_points()

  expect_equal(unname(unlist(pts[pts$i == 1 & pts$j == 0, c("x", "y", "z")])), c(5, 5, 5) + di)
  expect_equal(unname(unlist(pts[pts$i == 0 & pts$j == 1, c("x", "y", "z")])), c(5, 5, 5) + dj)
})

test_that("get_sample_points() is memoized and invalidated by setters", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(2, 2))
  pts1 <- s$get_sample_points()
  pts2 <- s$get_sample_points()
  expect_identical(pts1, pts2)

  s$set_origin(c(10, 0, 0))
  pts3 <- s$get_sample_points()
  expect_false(isTRUE(all.equal(pts1$x, pts3$x)))
  expect_equal(min(pts3$x), 10)
})

test_that("each setter invalidates the cache", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(2, 2))

  s$get_sample_points()
  s$set_direction(c(0, 1, 0), c(-1, 0, 0))
  expect_equal(unname(s$get_direction_i()), c(0, 1, 0))

  s$get_sample_points()
  s$set_spacing(c(5, 5))
  expect_equal(max(s$get_sample_points()$x %% 5), 0)

  s$get_sample_points()
  s$set_size(c(4, 4))
  expect_equal(nrow(s$get_sample_points()), 16)
})

test_that("translate() returns a new, independent SliceGeometry", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(2, 2))
  normal <- s$get_normal()
  s2 <- s$translate(5 * normal)

  expect_equal(s$get_origin(), c(0, 0, 0))
  expect_equal(s2$get_origin(), c(0, 0, 5))
  expect_equal(s2$get_direction(), s$get_direction())
  expect_equal(s2$get_spacing(), s$get_spacing())
  expect_equal(s2$get_size(), s$get_size())

  expect_error(s$translate(c(1, 2)), "length-3")
})

test_that("SliceGeometry$from_corners() derives direction_j, spacing, and size correctly", {
  p0 <- c(0, 0, 0)
  p1 <- c(3, 4, 0)
  di <- c(1, 0, 0)

  s <- SliceGeometry$from_corners(p0, p1, di, size = c(4, 5))
  expect_equal(unname(s$get_direction_j()), c(0, 1, 0))
  expect_equal(s$get_extent(), c(3, 4))
  expect_equal(s$get_size(), c(4L, 5L))

  pts <- s$get_sample_points()
  last <- pts[pts$i == max(pts$i) & pts$j == max(pts$j), c("x", "y", "z")]
  expect_equal(unname(unlist(last)), p1)

  s_spacing <- SliceGeometry$from_corners(p0, p1, di, spacing = c(1, 1))
  expect_equal(s_spacing$get_size(), c(4L, 5L))
})

test_that("SliceGeometry$from_corners() validates its inputs", {
  p0 <- c(0, 0, 0)
  p1 <- c(3, 4, 0)
  di <- c(1, 0, 0)

  expect_error(SliceGeometry$from_corners(p0, p1, di, spacing = c(1, 1), size = c(4, 5)), "exactly one")
  expect_error(SliceGeometry$from_corners(p0, p1, di), "exactly one")
  expect_error(SliceGeometry$from_corners(p0, p1, c(-1, 0, 0), size = c(4, 5)), "must point from")
  expect_error(SliceGeometry$from_corners(p0, c(6, 0, 0), di, size = c(4, 5)), "collinear")
  expect_error(SliceGeometry$from_corners(p0, p1, di, size = c(1, 5)), ">= 2")
})

test_that("as_sitk_reference_image() matches the manual sample-point geometry", {
  skip_if_not_installed("SimpleITK")

  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 2), c(3, 4))
  ref <- s$as_sitk_reference_image(spacing_k = 0.5)

  expect_equal(ref$GetSize(), c(3L, 4L, 1L))
  expect_equal(ref$GetOrigin(), c(0, 0, 0))
  expect_equal(ref$GetSpacing(), c(1, 2, 0.5))

  pt <- ref$TransformIndexToPhysicalPoint(c(2L, 3L, 0L))
  expected <- unname(unlist(s$get_sample_points()[s$get_sample_points()$i == 2 & s$get_sample_points()$j == 3, c("x", "y", "z")]))
  expect_equal(pt, expected)
})

test_that("print() runs without error and is invisible", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(2, 2))
  expect_invisible(s$print())
  expect_output(print(s), "<SliceGeometry>")
  expect_output(print(s), "normal:")
  expect_output(print(s), "extent:")
})

test_that("set_center() repositions the rectangle, keeping direction/spacing/size", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))
  s$set_center(c(10, 10, 0))

  expect_equal(unname(s$get_center()), c(10, 10, 0))
  expect_equal(s$get_origin(), c(9, 9, 0))
  expect_equal(s$get_spacing(), c(1, 1))
  expect_equal(s$get_size(), c(3L, 3L))

  expect_error(s$set_center(c(1, 2)), "length-3")
})

test_that("set_spacing_i()/set_spacing_j() update one axis independently", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))
  s$set_spacing_i(2)
  expect_equal(s$get_spacing(), c(2, 1))
  s$set_spacing_j(3)
  expect_equal(s$get_spacing(), c(2, 3))

  expect_error(s$set_spacing_i(0), "> 0")
  expect_error(s$set_spacing_j(-1), "> 0")
})

test_that("set_size_i()/set_size_j() update one axis independently", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))
  s$set_size_i(5)
  expect_equal(s$get_size(), c(5L, 3L))
  s$set_size_j(6)
  expect_equal(s$get_size(), c(5L, 6L))

  expect_error(s$set_size_i(0), ">= 1")
  expect_error(s$set_size_j(2.5), ">= 1")
})

test_that("set_extent() derives spacing or size as requested", {
  s <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(5, 5))

  s$set_extent(c(8, 8), adjust = "spacing")
  expect_equal(s$get_spacing(), c(2, 2))
  expect_equal(s$get_size(), c(5L, 5L))

  s$set_extent(c(8, 8), adjust = "size")
  expect_equal(s$get_size(), c(5L, 5L))

  expect_error(s$set_extent(c(-1, 1)), "non-negative")

  s_single <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(1, 3))
  expect_error(s_single$set_extent(c(4, 4), adjust = "spacing"), "size.*< 2 along an axis|Cannot derive")
})

test_that("SliceGeometry$from_center() builds a slice centered on the given point", {
  s <- SliceGeometry$from_center(c(1, 2, 3), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))
  expect_equal(unname(s$get_center()), c(1, 2, 3))
  expect_equal(s$get_size(), c(3L, 3L))
})

test_that("SliceGeometry$from_normal() derives an orthonormal basis with the requested normal", {
  s <- SliceGeometry$from_normal(origin = c(0, 0, 0), normal = c(0, 0, 5), spacing = c(1, 1), size = c(3, 3))
  expect_equal(unname(s$get_normal()), c(0, 0, 1))

  # explicit seed direction that is NOT already orthogonal to the normal
  s2 <- SliceGeometry$from_normal(
    origin = c(0, 0, 0), normal = c(1, 0, 0), spacing = c(1, 1), size = c(3, 3),
    direction_i = c(1, 1, 0)
  )
  expect_equal(unname(s2$get_normal()), c(1, 0, 0))
  expect_equal(sum(s2$get_direction_i() * s2$get_normal()), 0)

  expect_error(SliceGeometry$from_normal(c(0, 0, 0), c(0, 0, 0), c(1, 1), c(3, 3)), "nonzero")
  expect_error(
    SliceGeometry$from_normal(c(0, 0, 0), c(0, 0, 1), c(1, 1), c(3, 3), direction_i = c(0, 0, 2)),
    "parallel"
  )
})

test_that("SliceGeometry$from_image_axis() matches the source image's own geometry", {
  skip_if_not_installed("SimpleITK")

  img <- SimpleITK::Image(10L, 10L, 10L, "sitkFloat32")
  img$SetOrigin(c(-5, -5, -5))
  img$SetSpacing(c(1, 1, 1))

  s <- SliceGeometry$from_image_axis(img, "z", 0)
  expect_equal(s$get_size(), c(10L, 10L))
  expect_equal(s$get_spacing(), c(1, 1))
  expect_equal(s$get_origin(), c(-5, -5, 0))

  s_alias <- SliceGeometry$from_image_axis(img, "axial", 3.4)
  expect_equal(s_alias$get_origin()[3], 3)

  expect_error(SliceGeometry$from_image_axis(img, "bogus", 0), "Invalid `axis`")

  img2d <- SimpleITK::Image(c(4L, 4L), "sitkFloat32")
  expect_error(SliceGeometry$from_image_axis(img2d, "z", 0), "3D image")
})

test_that("SliceGeometry$from_image_axis() accepts numeric axis input too", {
  skip_if_not_installed("SimpleITK")

  img <- SimpleITK::Image(10L, 10L, 10L, "sitkFloat32")
  img$SetOrigin(c(-5, -5, -5))
  img$SetSpacing(c(1, 1, 1))

  s_word <- SliceGeometry$from_image_axis(img, "axial", 0)
  s_num <- SliceGeometry$from_image_axis(img, 3, 0)
  expect_equal(s_word$get_origin(), s_num$get_origin())
  expect_equal(s_word$get_size(), s_num$get_size())
})

test_that("SliceGeometry$from_bounds() round-trips get_bounds() and validates consistency", {
  s <- SliceGeometry$new(c(1, 2, 3), c(1, 0, 0), c(0, 1, 0), c(1, 2), c(4, 5))
  b <- s$get_bounds()

  s2 <- SliceGeometry$from_bounds(b, size = s$get_size())
  expect_equal(s2$get_origin(), s$get_origin())
  expect_equal(unname(s2$get_direction_i()), unname(s$get_direction_i()))
  expect_equal(unname(s2$get_direction_j()), unname(s$get_direction_j()))
  expect_equal(s2$get_size(), s$get_size())

  s3 <- SliceGeometry$from_bounds(b, spacing = s$get_spacing())
  expect_equal(s3$get_size(), s$get_size())

  expect_error(SliceGeometry$from_bounds(b, spacing = c(1, 1), size = c(4, 5)), "exactly one")
  expect_error(SliceGeometry$from_bounds(b[b$corner != "i1_j0", ], size = c(4, 5)), "must contain rows")

  b_bad <- b
  b_bad$x[b_bad$corner == "i1_j1"] <- 999
  expect_error(SliceGeometry$from_bounds(b_bad, size = s$get_size()), "inconsistent")

  b_degenerate <- b
  b_degenerate[b_degenerate$corner == "i1_j0", c("x", "y", "z")] <- b_degenerate[b_degenerate$corner == "i0_j0", c("x", "y", "z")]
  expect_error(SliceGeometry$from_bounds(b_degenerate, size = c(4, 5)), "Degenerate")
})
