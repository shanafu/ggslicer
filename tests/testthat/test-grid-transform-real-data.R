test_that("a real ANTs warp visibly distorts a grid built on the fixed image (mouse_4)", {
  skip_if_not_installed("SimpleITK")
  skip_if_no_testdata()

  base <- file.path(testdata_dir(), "mouse_4")
  fixed <- ReadImage_fix(file.path(base, "fmri_template_on_ccfv3_200um.nii.gz"))
  size <- fixed$GetSize()
  center_idx <- as.integer(c(size[1] %/% 2, size[2] %/% 2, size[3] %/% 2))
  center_world <- fixed$TransformIndexToPhysicalPoint(center_idx)

  geom <- SliceGeometry$from_image_axis(fixed, "axial", center_world[3])
  grid_df <- slice_grid(geometry = geom)

  # Capture what's needed from the warp image before constructing the
  # DisplacementFieldTransform -- doing so moves/invalidates the source Image
  # object (confirmed directly; see CLAUDE.md), though nothing here needs it
  # afterward anyway.
  inv_warp_img <- SimpleITK::ReadImage(
    file.path(base, "fmri_template_to_ccfv31InverseWarp.nii.gz"), "sitkVectorFloat64"
  )
  inv_transform <- SimpleITK::DisplacementFieldTransform(inv_warp_img)

  warped <- transform_points(grid_df, inv_transform)

  # A "grid" line with grid_axis == "i" has, by construction, zero variance in
  # its projection onto direction_i before warping (that's what "fixed local i"
  # means). After a genuine nonlinear warp, that projection should no longer be
  # constant for at least several lines (edge lines near/outside the brain may
  # legitimately see near-zero local deformation).
  origin <- geom$get_origin()
  di <- geom$get_direction_i()

  is_grid_i <- grid_df$part == "grid" & grid_df$grid_axis == "i"
  grid_i <- grid_df[is_grid_i, ]
  warped_i <- warped[is_grid_i, ] # transform_points() preserves row order/length exactly

  line_ids <- unique(grid_i$line_id)
  proj_std_after <- vapply(line_ids, function(lid) {
    sel <- grid_i$line_id == lid
    proj <- as.matrix(warped_i[sel, c("x", "y", "z")]) %*% di - sum(origin * di)
    stats::sd(proj)
  }, numeric(1))

  expect_true(max(proj_std_after) > 0.01)
})
