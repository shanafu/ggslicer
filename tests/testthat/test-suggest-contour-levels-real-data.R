test_that("suggest_contour_levels() gives sane, sorted, in-range levels on mouse_1's real average template", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  average <- ReadImage_fix(file.path(base, "DSURQE_40micron_average.mnc"))
  arr <- as.vector(SimpleITK::as.array(average))
  rng <- range(arr)

  quantile_levels <- suggest_contour_levels(average, method = "quantile", n = 5)
  expect_length(quantile_levels, 5)
  expect_equal(quantile_levels, sort(quantile_levels))
  expect_true(all(quantile_levels >= rng[1] & quantile_levels <= rng[2]))

  trough_levels <- suggest_contour_levels(average, method = "troughs", min_n = 1)
  expect_gte(length(trough_levels), 1)
  expect_equal(trough_levels, sort(trough_levels))
  expect_true(all(trough_levels >= rng[1] & trough_levels <= rng[2]))
})

test_that("suggest_contour_levels() rejects mouse_1's labels and mask (image and column form)", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  average <- ReadImage_fix(file.path(base, "DSURQE_40micron_average.mnc"))
  mask <- ReadImage_fix(file.path(base, "DSURQE_40micron_mask.mnc"))
  labels <- ReadImage_fix(file.path(base, "DSURQE_40micron_labels.mnc"))

  # image form: mask is caught by the data-driven check (integer-valued, only
  # 2 unique values). labels is NOT caught in image form -- a raw image has no
  # name to check, and DSURQE's atlas has 337 distinct integer labels, above
  # .looks_discrete()'s max_unique=50 heuristic threshold. This is a confirmed,
  # real limitation of the data-driven-only fallback (see @details), not a bug:
  # the name-based check (below, via a data frame column) is what actually
  # catches a labels image in practice.
  expect_error(suggest_contour_levels(mask, method = "quantile"), "discrete")
  expect_no_error(suggest_contour_levels(labels, method = "quantile"))

  # column form: name-based check fires first, before any data is even inspected.
  out <- slice_image(
    average, axis = "coronal", coordinate = 0,
    extra_images = list(mask = mask, labels = labels)
  )
  expect_error(suggest_contour_levels(out, column = "mask"), "discrete")
  expect_error(suggest_contour_levels(out, column = "labels"), "discrete")
  expect_no_error(suggest_contour_levels(out, column = "value", method = "quantile", n = 3))
})
