test_that("slice_label_contours() traces a real structure's boundary from human_1's annotation", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "human_1")
  annotation <- ReadImage_fix(file.path(base, "annotation.nii.gz"))

  # Label 10557 is the largest nonzero label present at z=8.0 (confirmed by direct
  # inspection: 13071 voxels in that slice, next-largest is ~1500).
  out <- slice_label_contours(annotation, axis = "axial", coordinate = 8.0, labels = 10557)

  expect_true(nrow(out) > 0)
  expect_setequal(unique(out$label), 10557)
  expect_true(all(c("package", "k", "label", "obj", "vertex", "x", "y", "z") %in% names(out)))
})

test_that("mask_fill = 'nan' vs 'zero' produce a real, dramatic difference on mouse_1", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  average <- ReadImage_fix(file.path(base, "DSURQE_40micron_average.mnc"))
  mask <- ReadImage_fix(file.path(base, "DSURQE_40micron_mask.mnc"))

  # At a very low level (0.5), real in-mask brain tissue intensity never crosses it
  # (confirmed by direct inspection), so any contour at this level under mask_fill =
  # "zero" is purely a fabricated edge at the mask boundary; mask_fill = "nan" should
  # then find no contour at all at this level.
  out_zero <- slice_contours(average, axis = "coronal", coordinate = 0, levels = 0.5, mask = mask, mask_fill = "zero")
  out_nan <- slice_contours(average, axis = "coronal", coordinate = 0, levels = 0.5, mask = mask, mask_fill = "nan")

  expect_true(nrow(out_zero) > 0)
  expect_equal(nrow(out_nan), 0)
})

test_that("slice_contours() finds a real iso-intensity boundary inside mouse_1's brain mask", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  average <- ReadImage_fix(file.path(base, "DSURQE_40micron_average.mnc"))
  mask <- ReadImage_fix(file.path(base, "DSURQE_40micron_mask.mnc"))

  # A mid-range level (500, well within the average's confirmed 0-2970 intensity
  # range) should produce real contours under both mask_fill modes.
  out_zero <- slice_contours(average, axis = "coronal", coordinate = 0, levels = 500, mask = mask, mask_fill = "zero")
  out_nan <- slice_contours(average, axis = "coronal", coordinate = 0, levels = 500, mask = mask, mask_fill = "nan")

  expect_true(nrow(out_zero) > 0)
  expect_true(nrow(out_nan) > 0)
})
