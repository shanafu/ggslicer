# Regression test for a specific question: does applying a MINC (.xfm-derived)
# transform to world points coming from ReadImage_fix()-corrected images suffer
# from the same x/y axis-flip issue that orientation_correction() corrects for
# when *reading* MINC images?
#
# Investigated directly (not assumed) using mouse_5's real registration outputs,
# which conveniently include ground truth: DSURQE_40micron_labels_on_CCFv3_25um.mnc
# is mouse_1's DSURQE labels already correctly warped onto CCFv3 space by the
# real MINC registration pipeline. Finding: NO axis-flip bug -- ITK's physical-
# point machinery (origin/direction/spacing) is self-consistent, so a world
# coordinate means the same physical location whether it came from a raw or a
# ReadImage_fix()-corrected reading of a MINC file; a transform defined in terms
# of physical points (an affine matrix, or a DisplacementFieldTransform built
# from a *plainly* read displacement-field volume -- never pass it through
# ReadImage_fix()) works correctly on either. This test guards that property
# with real data: several distinct anatomical labels from mouse_1, each
# transformed through mouse_5's real affine+grid transform, must land on the
# *same* label in the real, independently-computed ground-truth CCFv3-space
# labels.
test_that("a MINC transform applied to ReadImage_fix() points has no axis-flip bug", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  m1 <- file.path(testdata_dir(), "mouse_1")
  m5 <- file.path(testdata_dir(), "mouse_5")

  # Everything here is read exactly the way ggslicer users would read it.
  src_labels <- ReadImage_fix(file.path(m1, "DSURQE_40micron_labels.mnc"))
  target_labels <- ReadImage_fix(file.path(m5, "DSURQE_40micron_labels_on_CCFv3_25um.mnc"))

  # A pure single-Linear-block .xfm parses correctly via SimpleITK::ReadTransform()
  # (confirmed separately -- it's specifically multi-block/Grid_Transform .xfm
  # files that ReadTransform() silently mis-parses).
  linear <- SimpleITK::ReadTransform(file.path(m5, "affine", "MICe_DSURQE_affine.xfm"))

  # The grid (displacement field) volume must be read *plainly*, never through
  # ReadImage_fix() -- this is the one part of the pipeline that would actually
  # break if handled naively (see CLAUDE.md for why, and why it turns out not
  # to matter for the *lookup* itself but does matter as a documented rule).
  grid_img <- SimpleITK::ReadImage(file.path(m5, "MICe_DSURQE_grid_0.mnc"), "sitkVectorFloat64")
  grid <- SimpleITK::DisplacementFieldTransform(grid_img)

  # File order is [Linear, Grid_Transform] (apply Linear first, then Grid);
  # SimpleITK::CompositeTransform (R has no list-based constructor, unlike
  # Python's CompositeTransform([...]) -- use AddTransform() instead) applies
  # the *last*-added transform first, so add them in reverse: grid, then linear.
  composite <- SimpleITK::CompositeTransform(3L)
  composite$AddTransform(grid)
  composite$AddTransform(linear)

  arr <- SimpleITK::as.array(src_labels) # (i, j, k) order, corrected convention
  tab <- table(arr)
  tab <- tab[names(tab) != "0"]
  ordered_labels <- as.numeric(names(sort(tab, decreasing = TRUE)))
  nonzero_labels <- ordered_labels[1:8]
  expect_length(nonzero_labels, 8) # sanity: fixture actually has this many distinct labels

  mismatches <- list()
  for (label_value in nonzero_labels) {
    coords <- which(arr == label_value, arr.ind = TRUE)
    mid <- coords[ceiling(nrow(coords) / 2), ]
    idx_xyz <- as.integer(mid - 1) # which() is 1-indexed; SimpleITK indices are 0-indexed

    world_point <- src_labels$TransformIndexToPhysicalPoint(idx_xyz)
    transformed <- composite$TransformPoint(world_point)

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
