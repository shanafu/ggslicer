test_that("read_minc_transform() on real mouse_5 fixtures lands real anatomy on the correct label", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  m1 <- file.path(testdata_dir(), "mouse_1")
  m5 <- file.path(testdata_dir(), "mouse_5")

  src_labels <- ReadImage_fix(file.path(m1, "DSURQE_40micron_labels.mnc"))
  target_labels <- ReadImage_fix(file.path(m5, "DSURQE_40micron_labels_on_CCFv3_25um.mnc"))

  # The multi-block, Grid_Transform-containing top-level file -- exactly the
  # kind SimpleITK::ReadTransform() cannot parse correctly (confirmed
  # separately; see CLAUDE.md).
  xfm <- read_minc_transform(file.path(m5, "MICe_DSURQE.xfm"))
  expect_true(inherits(xfm, "_p_itk__simple__CompositeTransform"))

  arr <- SimpleITK::as.array(src_labels)
  tab <- table(arr)
  tab <- tab[names(tab) != "0"]
  labels <- as.numeric(names(sort(tab, decreasing = TRUE)))[1:8]

  mismatches <- list()
  for (label_value in labels) {
    coords <- which(arr == label_value, arr.ind = TRUE)
    mid <- coords[ceiling(nrow(coords) / 2), ]
    idx_xyz <- as.integer(mid - 1)

    world_point <- src_labels$TransformIndexToPhysicalPoint(idx_xyz)
    transformed <- xfm$TransformPoint(world_point)

    got <- tryCatch(
      target_labels$GetPixel(target_labels$TransformPhysicalPointToIndex(transformed)),
      error = function(e) NA_real_
    )
    if (!isTRUE(all.equal(got, as.numeric(label_value)))) {
      mismatches[[length(mismatches) + 1]] <- list(label = label_value, got = got)
    }
  }

  expect_length(mismatches, 0)
})

test_that("read_minc_transform() on a real single-block affine .xfm matches SimpleITK::ReadTransform()", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  m5 <- file.path(testdata_dir(), "mouse_5")
  path <- file.path(m5, "affine", "MICe_DSURQE_affine.xfm")

  # corrected = FALSE: comparing the raw parser's output against
  # SimpleITK::ReadTransform()'s own (also raw, unconjugated) parsing.
  via_parser <- read_minc_transform(path, corrected = FALSE)
  via_read_transform <- SimpleITK::ReadTransform(path)

  expect_equal(via_parser$GetParameters(), via_read_transform$GetParameters(), tolerance = 1e-8)
})
