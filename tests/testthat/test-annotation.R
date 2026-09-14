two_blob_label_image <- function(n = 60) {
  filled_sitk_image(c(n, n, 3), function(idx) {
    if ((idx[1] - 15)^2 + (idx[2] - 15)^2 <= 64) return(1)
    if ((idx[1] - 45)^2 + (idx[2] - 15)^2 <= 64) return(1)
    if ((idx[1] - 30)^2 + (idx[2] - 45)^2 <= 100) return(2)
    0
  })
}

test_that("contour_centroids() finds one centroid per connected component, not per label", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  contours <- slice_label_contours(img, axis = "z", coordinate = 0)

  out <- contour_centroids(contours)
  expect_true(all(c("package", "k", "label", "name", "obj", "x", "y", "z", "area") %in% names(out)))
  expect_equal(nrow(out), 3) # label 1 has 2 components, label 2 has 1

  label1 <- out[out$label == 1, ]
  expect_equal(nrow(label1), 2)
  expect_setequal(round(label1$x), c(15, 45))
  expect_true(all(abs(label1$y - 15) < 1e-6))

  label2 <- out[out$label == 2, ]
  expect_equal(nrow(label2), 1)
  expect_equal(round(label2$x), 30)
  expect_equal(round(label2$y), 45)
})

test_that("computed areas approximate the true circle areas", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  contours <- slice_label_contours(img, axis = "z", coordinate = 0)
  out <- contour_centroids(contours)

  # radius^2 = 64 -> area ~ pi*64 ~ 201; radius^2 = 100 -> area ~ pi*100 ~ 314
  expect_true(all(abs(out$area[out$label == 1] - pi * 64) < 15))
  expect_true(all(abs(out$area[out$label == 2] - pi * 100) < 15))
})

test_that("label_names maps label values to display text, falling back to the numeric id when unmapped", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  contours <- slice_label_contours(img, axis = "z", coordinate = 0)

  named <- contour_centroids(contours, label_names = c("1" = "RegionA", "2" = "RegionB"))
  expect_true(all(named$name[named$label == 1] == "RegionA"))
  expect_true(all(named$name[named$label == 2] == "RegionB"))

  unnamed <- contour_centroids(contours)
  expect_true(all(unnamed$name[unnamed$label == 1] == "1"))
  expect_true(all(unnamed$name[unnamed$label == 2] == "2"))

  partial <- contour_centroids(contours, label_names = c("1" = "RegionA"))
  expect_true(all(partial$name[partial$label == 1] == "RegionA"))
  expect_true(all(partial$name[partial$label == 2] == "2"))
})

test_that("min_area filters out small connected components", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  contours <- slice_label_contours(img, axis = "z", coordinate = 0)

  all_out <- contour_centroids(contours)
  expect_equal(nrow(all_out), 3)

  # label 1's components have area ~201; label 2's ~314 -- a threshold between
  # them should keep only label 2's single component.
  filtered <- contour_centroids(contours, min_area = 250)
  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$label, 2)

  # a threshold above everything drops all rows but keeps the right columns/types
  none <- contour_centroids(contours, min_area = 100000)
  expect_equal(nrow(none), 0)
  expect_true(all(c("package", "k", "label", "name", "obj", "x", "y", "z", "area") %in% names(none)))
})

test_that("contour_centroids() validates its arguments", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  contours <- slice_label_contours(img, axis = "z", coordinate = 0)

  expect_error(contour_centroids(data.frame(x = 1, y = 1)), "columns")
  expect_error(contour_centroids(contours, min_area = -1))
})

test_that("slice_label_annotations() matches calling slice_label_contours() + contour_centroids() manually", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()

  manual <- contour_centroids(
    slice_label_contours(img, axis = "z", coordinate = 0),
    label_names = c("1" = "RegionA", "2" = "RegionB"), min_area = 10
  )
  direct <- slice_label_annotations(
    img, axis = "z", coordinate = 0,
    label_names = c("1" = "RegionA", "2" = "RegionB"), min_area = 10
  )
  expect_equal(manual, direct)
})

test_that("slice_label_annotations() passes min_vertices through to slice_label_contours()", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  out <- slice_label_annotations(img, axis = "z", coordinate = 0, min_vertices = 100000)
  expect_equal(nrow(out), 0)
})

test_that("a thin arc ('C'-shaped) region's true-area centroid lands within its own band, not at the arc's center of curvature", {
  skip_if_not_installed("SimpleITK")
  # A thin ~80-degree arc slice of an annulus (radius 15-25) -- unlike a
  # nearly-*complete* ring (whose true area centroid is mathematically at
  # the ring's own geometric center, precisely because a symmetric ring's
  # mass is evenly distributed all the way around it -- confirmed directly:
  # an earlier version of this test used a ring missing only a small ~34-
  # degree wedge and got a centroid just ~2 units from center, which is
  # correct, not a bug, since that shape is nearly a full ring), a genuine
  # one-sided arc has no such symmetry, so its centroid should land
  # somewhere within its own 15-25 radius band.
  n <- 60
  cx <- 30; cy <- 30
  img <- filled_sitk_image(c(n, n, 3), function(idx) {
    dx <- idx[1] - cx; dy <- idx[2] - cy
    r <- sqrt(dx^2 + dy^2)
    ang <- atan2(dy, dx)
    if (r >= 15 && r <= 25 && ang > -0.7 && ang < 0.7) 1 else 0
  })
  contours <- slice_label_contours(img, axis = "z", coordinate = 0)
  out <- contour_centroids(contours)
  expect_equal(nrow(out), 1)

  dist_from_center <- sqrt((out$x - cx)^2 + (out$y - cy)^2)
  expect_true(dist_from_center > 10 && dist_from_center < 25)
})

test_that("slice_label_annotations_layer() builds a usable ggplot2 layer", {
  skip_if_not_installed("SimpleITK")
  img <- two_blob_label_image()
  out <- slice_label_annotations(img, axis = "z", coordinate = 0)

  layer <- slice_label_annotations_layer(out)
  expect_s3_class(layer, "ggproto")

  p <- ggplot2::ggplot(out, ggplot2::aes(x = x, y = y)) + layer
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})
