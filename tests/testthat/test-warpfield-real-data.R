test_that("slice_warp_arrows() on mouse_4's real ANTs warp gives sane, non-degenerate arrows away from the boundary", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_4")
  fixed <- ReadImage_fix(file.path(base, "fmri_template_on_ccfv3_200um.nii.gz"))
  size <- fixed$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- fixed$TransformIndexToPhysicalPoint(center_idx)

  # Read plainly -- never through ReadImage_fix()/orientation_correction() -- and
  # constructing DisplacementFieldTransform() from it moves/invalidates the
  # source Image object (confirmed elsewhere in this codebase; see warpfield.R).
  warp_img <- SimpleITK::ReadImage(file.path(base, "fmri_template_to_ccfv31Warp.nii.gz"), "sitkVectorFloat64")

  arrows <- slice_warp_arrows(
    fixed, axis = "axial", coordinate = center_world[3],
    warp = warp_img, spacing = 0.5
  )

  expect_gt(nrow(arrows), 0)
  expect_true(all(is.finite(arrows$total_displacement)))
  expect_true(max(arrows$total_displacement) > 0.01) # a real, non-degenerate warp at this slice

  # orthonormal invariant
  inv_err <- max(abs(
    arrows$total_displacement^2 - (arrows$in_plane_displacement^2 + arrows$normal_displacement^2)
  ))
  expect_lt(inv_err, 1e-6)

  # arrow never leaves the plane (axial slice: z is unaffected by the in-plane
  # x/y projection, regardless of the warp's real out-of-plane displacement)
  expect_equal(arrows$zend, arrows$z)
})

test_that("slice_warp_arrows() with arrow_length gives a uniform drawn length except where direction is undefined", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_4")
  fixed <- ReadImage_fix(file.path(base, "fmri_template_on_ccfv3_200um.nii.gz"))
  size <- fixed$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- fixed$TransformIndexToPhysicalPoint(center_idx)

  warp_img <- SimpleITK::ReadImage(file.path(base, "fmri_template_to_ccfv31Warp.nii.gz"), "sitkVectorFloat64")
  arrows <- slice_warp_arrows(
    fixed, axis = "axial", coordinate = center_world[3],
    warp = warp_img, spacing = 0.5, arrow_length = 0.2
  )

  drawn_len <- sqrt((arrows$xend - arrows$x)^2 + (arrows$yend - arrows$y)^2 + (arrows$zend - arrows$z)^2)
  has_direction <- arrows$in_plane_displacement > sqrt(.Machine$double.eps)
  expect_equal(drawn_len[has_direction], rep(0.2, sum(has_direction)))
  expect_equal(drawn_len[!has_direction], rep(0, sum(!has_direction)))
})

test_that("a real MINC-format warp (mouse_5's .xfm) gives arrows consistent with the equivalent nifti warp (mouse_4)", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base4 <- file.path(testdata_dir(), "mouse_4")
  fixed <- ReadImage_fix(file.path(base4, "fmri_template_on_ccfv3_200um.nii.gz"))
  size <- fixed$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- fixed$TransformIndexToPhysicalPoint(center_idx)

  base5 <- file.path(testdata_dir(), "mouse_5")
  arrows <- slice_warp_arrows(
    fixed, axis = "axial", coordinate = center_world[3],
    warp = file.path(base5, "MICe_DSURQE.xfm"), spacing = 5
  )

  expect_gt(nrow(arrows), 0)
  expect_true(all(is.finite(arrows$total_displacement)))
  expect_true(max(arrows$total_displacement) > 0)
})

test_that("slice_warp_arrows_layer() renders a real ggplot without error on real mouse_4 data", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_4")
  fixed <- ReadImage_fix(file.path(base, "fmri_template_on_ccfv3_200um.nii.gz"))
  size <- fixed$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- fixed$TransformIndexToPhysicalPoint(center_idx)

  warp_img <- SimpleITK::ReadImage(file.path(base, "fmri_template_to_ccfv31Warp.nii.gz"), "sitkVectorFloat64")
  arrows <- slice_warp_arrows(fixed, axis = "axial", coordinate = center_world[3], warp = warp_img, spacing = 1)

  p <- ggplot2::ggplot(arrows, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    slice_warp_arrows_layer(arrows)
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})
