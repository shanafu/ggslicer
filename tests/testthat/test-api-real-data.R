test_that("slice_image() with human_1's annotation as an extra image uses nearest-neighbor by name", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "human_1")
  template <- ReadImage_fix(file.path(base, "mni_icbm152_t1_tal_nlin_sym_09b_hires.nii"))
  # Cast to float: the raw annotation is uint32, and SimpleITK::Resample() preserves
  # the input pixel type, which would silently round a forced-linear sample back to
  # an integer and defeat the point of this test (checking that linear vs.
  # nearest-neighbor actually differ).
  annotation <- SimpleITK::Cast(ReadImage_fix(file.path(base, "annotation.nii.gz")), "sitkFloat32")

  # build_slice_geometry()/SliceGeometry$from_image_axis() snaps the in-plane sample
  # grid exactly onto the source image's own voxel grid, so on-grid sampling can never
  # distinguish nearest-neighbor from linear (both land exactly on a voxel center).
  # Shift the in-plane origin by a quarter voxel (via the geometry= escape hatch) so
  # samples genuinely fall between annotation's voxel centers, crossing the many label
  # boundaries confirmed (by direct inspection) to exist at z=8.0.
  on_grid <- build_slice_geometry(template, "axial", 8.0)
  base_slice <- on_grid$get_packages()$package_1$get_base_slice()
  quarter_voxel <- 0.25 * (
    base_slice$get_direction_i() * base_slice$get_spacing()[1] +
      base_slice$get_direction_j() * base_slice$get_spacing()[2]
  )
  off_grid <- base_slice$translate(quarter_voxel)

  out <- slice_image(template, geometry = off_grid, extra_images = list(annotation = annotation))
  vals <- out$annotation[!is.na(out$annotation)]
  expect_true(all(abs(vals - round(vals)) < 1e-6))

  out_linear <- slice_image(
    template, geometry = off_grid,
    extra_images = list(annotation = annotation),
    interpolator_overrides = list(annotation = "sitkLinear")
  )
  vals_linear <- out_linear$annotation[!is.na(out_linear$annotation)]
  expect_true(any(abs(vals_linear - round(vals_linear)) > 1e-6))
})

test_that("slice_image() works correctly on both MINC- and NIfTI-loaded versions of the same template", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  # Note: ReadImage_fix()'s orientation_correction() is a display-orientation
  # fix (it mirrors the pixel data and derives a new origin from the volume's
  # own extent), not a true world-coordinate registration between formats. For
  # this fixture, the MINC-corrected origin and the NIfTI file's own stored
  # origin genuinely differ by more than a rounding error along y (the brain
  # isn't vertically centered in the volume), so the two are not expected to
  # land on identical world coordinates for the same nominal slice — only
  # verified here is that both load into a consistent RAS+ (positive x/y)
  # orientation and sample correctly through the new high-level API.
  base <- file.path(testdata_dir(), "human_1")
  mnc <- ReadImage_fix(file.path(base, "mni_icbm152_t1_tal_nlin_sym_09b_hires.mnc"))
  nii <- ReadImage_fix(file.path(base, "mni_icbm152_t1_tal_nlin_sym_09b_hires.nii"))

  expect_true(all(mnc$GetOrigin()[1:2] > 0))
  expect_true(all(nii$GetOrigin()[1:2] > 0))

  out_mnc <- slice_image(mnc, axis = "axial", coordinate = mnc$GetOrigin()[3])
  out_nii <- slice_image(nii, axis = "axial", coordinate = nii$GetOrigin()[3])

  expect_equal(nrow(out_mnc), 394 * 466)
  expect_equal(nrow(out_nii), 394 * 466)
  expect_false(all(is.na(out_mnc$value)))
  expect_false(all(is.na(out_nii$value)))
})

test_that("mouse_1's mask column is usable for tidy-side filtering (dplyr::filter())", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_1")
  average <- ReadImage_fix(file.path(base, "DSURQE_40micron_average.mnc"))
  mask <- ReadImage_fix(file.path(base, "DSURQE_40micron_mask.mnc"))
  labels <- ReadImage_fix(file.path(base, "DSURQE_40micron_labels.mnc"))

  out <- slice_image(
    average, axis = "coronal", coordinate = 0,
    extra_images = list(mask = mask, labels = labels)
  )
  expect_true(all(c("mask", "labels") %in% names(out)))

  filtered <- dplyr::filter(out, mask > 0)
  expect_lt(nrow(filtered), nrow(out))
})

test_that("build_slice_geometry() places human_3's uneven slice coordinates correctly", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "human_3")
  template <- ReadImage_fix(file.path(base, "mni_icbm152_t1_tal_nlin_sym_09b_hires.nii"))

  pset <- build_slice_geometry(template, "coronal", c(9.0, 8.5))
  expect_length(pset$get_package_names(), 2)

  pts <- pset$get_sample_points()
  y_by_package <- tapply(pts$y, pts$package, function(v) unique(round(v, 4)))
  expect_setequal(unlist(y_by_package), c(9.0, 8.5))
})

test_that("sample_images() combines mouse_2's multi-resolution template/labels/gene-expression", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_2")
  labels <- ReadImage_fix(file.path(base, "AMBA_relabeled_backsampled_50um.mnc"))
  template25 <- ReadImage_fix(file.path(base, "average_template_25.mnc"))
  bdnf <- ReadImage_fix(file.path(base, "Bdnf_79587720.mnc"))

  geom <- SliceGeometry$from_image_axis(labels, "coronal", labels$GetOrigin()[2])
  out <- sample_images(geom, list(labels = labels, template = template25, bdnf = bdnf))

  expect_true(all(c("labels", "template", "bdnf") %in% names(out)))
  expect_false(all(is.na(out$template)))
  expect_false(all(is.na(out$bdnf)))
})
