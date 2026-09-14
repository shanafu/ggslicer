# Normalize `warp` into a SimpleITK::Transform object. Accepts an
# already-built SimpleITK::Transform (any subclass -- Affine, Composite,
# DisplacementFieldTransform, etc. -- used directly), a SimpleITK::Image with
# 3 components per pixel (a raw displacement field, e.g. an ANTs
# *Warp.nii.gz -- confirmed directly this is how such files actually load
# via SimpleITK::ReadImage(), as a genuine vector-pixel image, *not* a
# literal 4D scalar array -- auto-wrapped via
# SimpleITK::DisplacementFieldTransform()), or a file path to either (".xfm"
# routed through read_minc_transform(), reusing its Grid_Transform/Linear/
# composite-block handling; any other path first tried as a vector image,
# falling back to SimpleITK::ReadTransform() for a genuine transform-format
# file such as .tfm/.mat/.h5).
.as_warp_transform <- function(warp) {
  if (is.character(warp)) {
    if (tolower(tools::file_ext(warp)) == "xfm") {
      return(read_minc_transform(warp))
    }
    img <- tryCatch(SimpleITK::ReadImage(warp), error = function(e) NULL)
    if (!is.null(img)) {
      if (img$GetNumberOfComponentsPerPixel() != 3) {
        stop(
          "`warp` (\"", warp, "\") is an image but does not have 3 components ",
          "per pixel, so it can't be a displacement field.",
          call. = FALSE
        )
      }
      return(SimpleITK::DisplacementFieldTransform(img))
    }
    return(SimpleITK::ReadTransform(warp))
  }
  if (inherits(warp, "_p_itk__simple__Image")) {
    if (warp$GetNumberOfComponentsPerPixel() != 3) {
      stop("`warp` must be a 3-component vector image (a displacement field).", call. = FALSE)
    }
    return(SimpleITK::DisplacementFieldTransform(warp))
  }
  if (inherits(warp, "_p_itk__simple__Transform")) {
    return(warp)
  }
  stop(
    "`warp` must be a SimpleITK Transform object, a 3-component vector Image ",
    "(a displacement field), or a file path to either.",
    call. = FALSE
  )
}

# Build the arrow overlay for a single SliceGeometry: a coarser lattice of
# anchor points (same origin/direction_i/direction_j as slice_geom, but
# spaced `spacing` apart in both directions -- mirroring slice_grid()'s own
# spacing convention exactly, so the two overlays share one "how coarse is
# reasonable" default), each displaced by `transform` and decomposed into
# in-plane (drawn) and normal (metadata-only) components.
.slice_warp_arrows_one <- function(slice_geom, transform, spacing, arrow_length) {
  extent <- slice_geom$get_extent()
  if (is.null(spacing)) spacing <- min(extent) / 10

  size <- round(extent / spacing) + 1

  coarse <- SliceGeometry$new(
    origin = slice_geom$get_origin(),
    direction_i = slice_geom$get_direction_i(),
    direction_j = slice_geom$get_direction_j(),
    spacing = c(spacing, spacing),
    size = size
  )
  points <- coarse$get_sample_points()

  di <- slice_geom$get_direction_i()
  dj <- slice_geom$get_direction_j()
  normal <- slice_geom$get_normal()

  starts <- cbind(points$x, points$y, points$z)
  # SimpleITK has no vectorized/batch transform-point API (same constraint
  # documented for transform_points()); progress is reported the same way,
  # via pbapply::pbapply() -- silent by default outside an interactive
  # session, and here the point count is deliberately small (a sparse arrow
  # lattice, not a fine grid), so the per-point binding-overhead cost
  # documented in CLAUDE.md's performance section is not a practical concern.
  transformed <- pbapply::pbapply(starts, 1, function(p) transform$TransformPoint(p))

  d <- t(transformed) - starts
  d_i <- drop(d %*% di)
  d_j <- drop(d %*% dj)
  d_n <- drop(d %*% normal)

  in_plane_mag <- sqrt(d_i^2 + d_j^2)
  total_mag <- sqrt(in_plane_mag^2 + d_n^2)

  if (is.null(arrow_length)) {
    draw_i <- d_i
    draw_j <- d_j
  } else {
    # Direction is undefined where in-plane displacement is ~0 (pure
    # out-of-plane movement, or literally zero displacement) -- leave those
    # as a zero-length arrow rather than rescaling an arbitrary direction;
    # normal_displacement still reports the true (possibly nonzero) movement
    # at that point for optional color/alpha encoding.
    has_direction <- in_plane_mag > sqrt(.Machine$double.eps)
    scale <- ifelse(has_direction, arrow_length / in_plane_mag, 0)
    draw_i <- d_i * scale
    draw_j <- d_j * scale
  }

  ends <- starts + outer(draw_i, di) + outer(draw_j, dj)

  tibble::tibble(
    i = points$i, j = points$j,
    x = points$x, y = points$y, z = points$z,
    xend = ends[, 1], yend = ends[, 2], zend = ends[, 3],
    in_plane_displacement = in_plane_mag,
    normal_displacement = d_n,
    total_displacement = total_mag
  )
}

#' Build a per-slice overlay of nonlinear-warp displacement arrows
#'
#' @description
#' An alternative to warping a grid (`slice_grid()` + [transform_points()])
#' for visualizing a nonlinear registration warp: at a lattice of points
#' across a slice, shows the local displacement as an arrow rather than the
#' aggregate distortion of a warped grid line. Ready for
#' `ggplot2::geom_segment(aes(x = x, y = y, xend = xend, yend = yend))`.
#'
#' `warp` is normalized to one `SimpleITK::Transform` internally (see
#' `.as_warp_transform()`), so this works uniformly whether you have a raw
#' displacement-field image (e.g. an ANTs `*Warp.nii.gz`), a `.xfm` file, or
#' an already-loaded transform object of any kind (affine, composite,
#' displacement field) -- the displacement at each point is just
#' `transform$TransformPoint(point) - point`, exactly the mechanism
#' [transform_points()] already uses.
#'
#' Since plotting an arrow at every sampled voxel would be an illegible,
#' overlapping mess, arrows are placed on a coarser lattice (`spacing` apart
#' in both directions -- assumed equal, matching `slice_grid()`'s own
#' convention and default), not at the source image's native resolution.
#'
#' Each 3D displacement vector is decomposed, using the slice's own
#' orthonormal basis, into an in-plane part (`direction_i`/`direction_j`
#' components -- this is what's actually drawn, so the arrow always stays
#' exactly within the slice plane) and a normal (out-of-plane) part, which
#' cannot be drawn as an in-plane arrow but is still reported
#' (`normal_displacement`) so it can be mapped to e.g. `color`/`alpha` if a
#' user wants to flag where the true 3D movement is mostly out-of-plane.
#'
#' @param image A `SimpleITK` image, or a file path (read internally via
#'   [ReadImage_fix()]). Required unless `geometry` is supplied.
#' @param axis Which image axis is out-of-plane (see [build_slice_geometry()]).
#'   Required unless `geometry` is supplied; must not be supplied together
#'   with `geometry`.
#' @param coordinate One or more world coordinates along `axis` (see
#'   [build_slice_geometry()]). Required unless `geometry` is supplied; must
#'   not be supplied together with `geometry`.
#' @param warp A `SimpleITK::Transform` object, a `SimpleITK` image with 3
#'   components per pixel (a raw displacement field; auto-wrapped via
#'   `SimpleITK::DisplacementFieldTransform()`), or a file path to either
#'   (`.xfm` routed through [read_minc_transform()]). If a raw image object
#'   is passed directly (not a path), be aware that constructing a
#'   `DisplacementFieldTransform()` from it moves/invalidates the source
#'   image object (a confirmed `SimpleITK` behavior, also noted for
#'   `read_minc_transform()`'s `Grid_Transform` handling) -- re-read the
#'   image if you need to use it again afterward.
#' @param spacing Distance between adjacent arrow anchor points, the same for
#'   both axes. Defaults (if `NULL`) to `min(extent) / 10`, matching
#'   `slice_grid()`'s own default, so an arrow overlay and a grid overlay
#'   default to the same visual density.
#' @param arrow_length If `NULL` (default), each arrow's drawn length is the
#'   true in-plane displacement at that point (world units). If given (a
#'   single positive number, world units), every arrow's in-plane component
#'   is instead rescaled to exactly this length (direction preserved) --
#'   useful for visually emphasizing direction when displacement magnitude
#'   varies a lot across the slice. The true magnitudes
#'   (`in_plane_displacement`/`normal_displacement`/`total_displacement`) are
#'   always reported regardless of this setting. Points with ~zero in-plane
#'   displacement (direction undefined) are left as a zero-length arrow even
#'   when `arrow_length` is set.
#' @param invert If `TRUE`, apply the inverse of `warp` (`GetInverse()`) --
#'   see [transform_points()] for the same caveat regarding
#'   `DisplacementFieldTransform`/composite transforms.
#' @param geometry A pre-built [SliceGeometry]/[SlicePackage]/[SlicePackageSet],
#'   used instead of `axis`/`coordinate`. Must not be supplied together with
#'   `axis`/`coordinate`.
#' @return A tibble with columns `package`, `k`, `i`, `j` (lattice indices,
#'   0-indexed, at the coarse `spacing` resolution -- not the source image's
#'   own voxel grid), `x`, `y`, `z` (arrow start, world coordinates), `xend`,
#'   `yend`, `zend` (arrow end, world coordinates, always exactly in-plane),
#'   `in_plane_displacement`, `normal_displacement` (signed), and
#'   `total_displacement` (the full, true 3D displacement magnitude,
#'   `sqrt(in_plane_displacement^2 + normal_displacement^2)`, regardless of
#'   `arrow_length`).
#' @examples
#' \dontrun{
#' fixed_image <- ReadImage_fix("fixed_template.nii")
#' arrows <- slice_warp_arrows(
#'   fixed_image, axis = "axial", coordinate = 0,
#'   warp = "registration.xfm"
#' )
#'
#' library(ggplot2)
#' ggplot(arrows, aes(x = x, y = y, xend = xend, yend = yend)) +
#'   geom_segment(arrow = grid::arrow(length = grid::unit(0.1, "inches")))
#' }
#' @export
slice_warp_arrows <- function(image = NULL, axis = NULL, coordinate = NULL, warp,
                               spacing = NULL, arrow_length = NULL, invert = FALSE,
                               geometry = NULL) {
  if (!is.null(spacing) && (!is.numeric(spacing) || length(spacing) != 1 || spacing <= 0)) {
    stop("`spacing` must be a single positive number (or NULL for the default).", call. = FALSE)
  }
  if (!is.null(arrow_length) && (!is.numeric(arrow_length) || length(arrow_length) != 1 || arrow_length <= 0)) {
    stop("`arrow_length` must be a single positive number (or NULL for true displacement).", call. = FALSE)
  }

  transform <- .as_warp_transform(warp)
  if (invert) {
    transform <- transform$GetInverse()
  }

  package_set <- .resolve_geometry_only(image, axis, coordinate, geometry)

  results <- list()
  for (pkg_name in package_set$get_package_names()) {
    pkg <- package_set$get_packages()[[pkg_name]]
    n_k <- pkg$get_size_k()
    for (k in 0:(n_k - 1)) {
      slice_geom <- pkg$get_slice(k)
      arrows <- .slice_warp_arrows_one(slice_geom, transform, spacing, arrow_length)
      arrows$package <- pkg_name
      arrows$k <- k
      results[[length(results) + 1]] <- arrows
    }
  }

  out <- dplyr::bind_rows(results)
  dplyr::select(
    out, "package", "k", "i", "j", "x", "y", "z", "xend", "yend", "zend",
    "in_plane_displacement", "normal_displacement", "total_displacement"
  )
}

#' Ready-made ggplot2 layer for a warp-arrow overlay produced by slice_warp_arrows()
#'
#' @description
#' Wraps `ggplot2::geom_segment()` with an arrowhead (`grid::arrow()`), so a
#' `slice_warp_arrows()` data frame can be added to a plot directly, without
#' the caller having to remember the `arrow =` argument -- the same
#' ready-made-layer convenience `slice_grid_layers()` provides for `slice_grid()`.
#'
#' @param arrows_df A data frame as returned by [slice_warp_arrows()].
#' @param color,linewidth,alpha Style for the arrow segments.
#' @param arrow_length_inches Length of the drawn arrowhead itself (inches,
#'   via `grid::arrow()`'s `length`) -- unrelated to `slice_warp_arrows()`'s
#'   own `arrow_length` argument, which controls the *shaft*'s world-space
#'   length, not the arrowhead's on-page size.
#' @param arrow_type `"open"` (default) or `"closed"`, passed to `grid::arrow()`.
#' @return A single `ggplot2::geom_segment()` layer, with an arrowhead. Add it
#'   (e.g. `plt + layer`) to a plot that already maps `x`/`y` via its own `aes()`.
#' @examples
#' \dontrun{
#' fixed_image <- ReadImage_fix("fixed_template.nii")
#' arrows <- slice_warp_arrows(
#'   fixed_image, axis = "axial", coordinate = 0,
#'   warp = "registration.xfm"
#' )
#' layer <- slice_warp_arrows_layer(arrows)
#'
#' library(ggplot2)
#' ggplot(arrows, aes(x = x, y = y, xend = xend, yend = yend)) + layer
#' }
#' @importFrom rlang .data
#' @export
slice_warp_arrows_layer <- function(arrows_df, color = "#2C7FB8", linewidth = 0.5, alpha = 1,
                                     arrow_length_inches = 0.08, arrow_type = c("open", "closed")) {
  arrow_type <- match.arg(arrow_type)
  # Use the unqualified `.data` pronoun (imported via @importFrom rlang .data
  # below), not `rlang::.data` -- see the identical note in grid.R's
  # slice_grid_layers(): ggplot2's tidy-eval data mask only special-cases the
  # bare `.data` symbol, and a bare column symbol (`x = x`) instead would
  # trigger an R CMD check "no visible binding for global variable" NOTE.
  ggplot2::geom_segment(
    data = arrows_df,
    mapping = ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    color = color, linewidth = linewidth, alpha = alpha,
    arrow = grid::arrow(length = grid::unit(arrow_length_inches, "inches"), type = arrow_type)
  )
}
