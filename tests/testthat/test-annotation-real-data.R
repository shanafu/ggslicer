test_that("slice_label_annotations() on mouse_1's real DSURQE atlas gives one centroid per real connected component", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  labels <- ReadImage_fix(file.path(base, "DSURQE_40micron_labels.mnc"))
  size <- labels$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- labels$TransformIndexToPhysicalPoint(center_idx)

  contours <- slice_label_contours(labels, axis = "coronal", coordinate = center_world[2])
  out <- contour_centroids(contours)

  expect_gt(nrow(out), 0)
  # every (label, obj) pair present in the contours should have exactly one
  # centroid row, and vice versa -- no components lost or duplicated
  expect_setequal(
    paste(contours$label, contours$obj),
    paste(out$label, out$obj)
  )
  expect_true(all(is.finite(out$x) & is.finite(out$y) & is.finite(out$z)))
  expect_true(all(out$area > 0))

  # unmapped labels fall back to their own numeric id as text
  expect_true(all(out$name == as.character(as.integer(round(out$label)))))
})

test_that("min_area meaningfully reduces real, cluttered atlas output while keeping the largest regions", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  labels <- ReadImage_fix(file.path(base, "DSURQE_40micron_labels.mnc"))
  size <- labels$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- labels$TransformIndexToPhysicalPoint(center_idx)

  all_out <- slice_label_annotations(labels, axis = "coronal", coordinate = center_world[2])
  filtered <- slice_label_annotations(labels, axis = "coronal", coordinate = center_world[2], min_area = 1)

  expect_lt(nrow(filtered), nrow(all_out))
  expect_true(all(filtered$area >= 1))

  # the largest real component should always survive a modest area threshold
  largest <- all_out[which.max(all_out$area), ]
  expect_true(any(filtered$label == largest$label & filtered$obj == largest$obj))
})

test_that("label_names correctly annotates a handful of real DSURQE label ids", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  labels <- ReadImage_fix(file.path(base, "DSURQE_40micron_labels.mnc"))
  size <- labels$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- labels$TransformIndexToPhysicalPoint(center_idx)

  contours <- slice_label_contours(labels, axis = "coronal", coordinate = center_world[2])
  present_labels <- unique(contours$label)
  skip_if(length(present_labels) < 2, "not enough distinct real labels at this slice to test naming")

  name_map <- setNames(paste0("Region_", present_labels[1:2]), as.character(present_labels[1:2]))
  out <- contour_centroids(contours, label_names = name_map)

  expect_true(all(out$name[out$label == present_labels[1]] == paste0("Region_", present_labels[1])))
  expect_true(all(out$name[out$label == present_labels[2]] == paste0("Region_", present_labels[2])))
})

test_that("slice_label_annotations_layer() renders a real ggplot with real atlas annotations overlaid on the template", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  average <- ReadImage_fix(file.path(base, "DSURQE_40micron_average.mnc"))
  labels <- ReadImage_fix(file.path(base, "DSURQE_40micron_labels.mnc"))
  size <- labels$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- labels$TransformIndexToPhysicalPoint(center_idx)

  img_df <- slice_image(average, axis = "coronal", coordinate = center_world[2])
  annotations <- slice_label_annotations(
    labels, axis = "coronal", coordinate = center_world[2], min_area = 1
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = img_df, ggplot2::aes(x = x, y = y, fill = value)) +
    slice_label_annotations_layer(annotations, color = "white")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})
