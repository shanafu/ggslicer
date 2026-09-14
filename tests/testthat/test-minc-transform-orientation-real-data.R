# Regression test for a real bug and its fix: does a MINC (.xfm-derived)
# transform correctly apply to world points coming from ReadImage_fix()-
# corrected images?
#
# History: an earlier version of this file concluded "no axis-flip bug"
# based on this same 8-real-label round trip passing under the *old*,
# now-removed FlipImageFilter-based orientation_correction(). That
# conclusion was incomplete: the old correction's images were themselves
# wrong (confirmed separately -- see test-io.R and CLAUDE.md -- against
# real MINC/NIfTI ground truth via mincheader/fslhd), but happened to
# combine with *raw* (unconverted) .xfm transforms in a way that still
# produced correct results for this specific fixture's geometry -- not a
# general guarantee.
#
# The real, verified picture (see CLAUDE.md for the full investigation):
# a .xfm file's own matrix/displacement values are defined in MINC's
# *native* coordinate convention (confirmed directly: they only reproduce
# real registered anatomy when applied to points from *plainly*-read MINC
# images, not ReadImage_fix()-corrected ones). ReadImage_fix()-corrected
# images, however, use a different (RAS/LPS-consistent) convention. These
# two conventions are related by a fixed, known operation (negate x/y,
# leave z), so read_minc_transform()'s default (corrected = TRUE) conjugates
# the parsed transform by that same operation, making it correct for points
# from ReadImage_fix()-corrected images -- which is what every function in
# this package actually produces.
test_that("read_minc_transform()'s default correctly transforms points from ReadImage_fix()-corrected images", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  m1 <- file.path(testdata_dir(), "mouse_1")
  m5 <- file.path(testdata_dir(), "mouse_5")

  # Everything here is read exactly the way ggslicer users would read it.
  src_labels <- ReadImage_fix(file.path(m1, "DSURQE_40micron_labels.mnc"))
  target_labels <- ReadImage_fix(file.path(m5, "DSURQE_40micron_labels_on_CCFv3_25um.mnc"))

  # The single concatenated file (Linear + Grid_Transform, in that order);
  # read_minc_transform()'s default corrected = TRUE conjugates the whole
  # thing at once.
  xfm <- read_minc_transform(file.path(m5, "MICe_DSURQE.xfm"))

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

test_that("the same real transform also works end-to-end on *plainly*-read (uncorrected) MINC images with corrected = FALSE", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  m1 <- file.path(testdata_dir(), "mouse_1")
  m5 <- file.path(testdata_dir(), "mouse_5")

  # Plain reads -- never through ReadImage_fix() -- paired with
  # corrected = FALSE, i.e. the transform exactly as written in the file.
  src_labels <- SimpleITK::ReadImage(file.path(m1, "DSURQE_40micron_labels.mnc"))
  target_labels <- SimpleITK::ReadImage(file.path(m5, "DSURQE_40micron_labels_on_CCFv3_25um.mnc"))
  xfm <- read_minc_transform(file.path(m5, "MICe_DSURQE.xfm"), corrected = FALSE)

  arr <- SimpleITK::as.array(src_labels)
  tab <- table(arr)
  tab <- tab[names(tab) != "0"]
  nonzero_labels <- as.numeric(names(sort(tab, decreasing = TRUE)))[1:8]

  mismatches <- list()
  for (label_value in nonzero_labels) {
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

test_that("mixing corrected and uncorrected conventions gives wrong results (confirms the conjugation is load-bearing, not a no-op)", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  m1 <- file.path(testdata_dir(), "mouse_1")
  m5 <- file.path(testdata_dir(), "mouse_5")

  # ReadImage_fix()-corrected images, but the RAW (unconjugated) transform --
  # a genuine coordinate-convention mismatch.
  src_labels <- ReadImage_fix(file.path(m1, "DSURQE_40micron_labels.mnc"))
  target_labels <- ReadImage_fix(file.path(m5, "DSURQE_40micron_labels_on_CCFv3_25um.mnc"))
  xfm_raw <- read_minc_transform(file.path(m5, "MICe_DSURQE.xfm"), corrected = FALSE)

  arr <- SimpleITK::as.array(src_labels)
  tab <- table(arr)
  tab <- tab[names(tab) != "0"]
  nonzero_labels <- as.numeric(names(sort(tab, decreasing = TRUE)))[1:8]

  matches <- 0
  for (label_value in nonzero_labels) {
    coords <- which(arr == label_value, arr.ind = TRUE)
    mid <- coords[ceiling(nrow(coords) / 2), ]
    idx_xyz <- as.integer(mid - 1)
    world_point <- src_labels$TransformIndexToPhysicalPoint(idx_xyz)
    transformed <- xfm_raw$TransformPoint(world_point)
    got <- tryCatch(
      target_labels$GetPixel(target_labels$TransformPhysicalPointToIndex(transformed)),
      error = function(e) NA_real_
    )
    if (isTRUE(all.equal(got, as.numeric(label_value)))) matches <- matches + 1
  }
  expect_equal(matches, 0)
})
