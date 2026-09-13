#' @importFrom rlang :=
#' @importFrom magrittr %>%
NULL

# Reshape one slice's sampled intensity column (in the i-fastest, then j,
# then k row order that SlicePackage$get_sample_points()/sample_intensity()
# always produce) into an (n_i x n_j) matrix suitable for contourLines().
.slice_intensity_matrix <- function(values, n_i, n_j) {
  matrix(values, nrow = n_i, ncol = n_j)
}

# Apply mask_fill to a sampled intensity vector given the corresponding
# sampled mask vector (already nearest-neighbor sampled onto the same grid).
.apply_mask_fill <- function(values, mask_values, mask_fill) {
  excluded <- !is.na(mask_values) & mask_values < 0.5
  if (mask_fill == "zero") {
    values[excluded] <- 0
  } else {
    values[excluded] <- NaN
  }
  values
}

# Map a contourLines()-style path (local in-plane coordinates, in the same
# physical units as the x/y sequences passed to contourLines()) to world
# coordinates via a SliceGeometry's origin/direction_i/direction_j.
.contour_path_to_world <- function(path_x, path_y, slice_geom) {
  origin <- slice_geom$get_origin()
  di <- slice_geom$get_direction_i()
  dj <- slice_geom$get_direction_j()
  n_pts <- length(path_x)
  world <- matrix(origin, nrow = n_pts, ncol = 3, byrow = TRUE) +
    outer(path_x, di) + outer(path_y, dj)
  world
}

# Validate/resolve the shared geometry + axis/coordinate + image/mask
# arguments common to slice_contours()/slice_label_contours(). Returns a
# list(image=, mask=, package_set=).
.resolve_contour_inputs <- function(image, axis, coordinate, mask, geometry) {
  has_geometry <- !is.null(geometry)
  if (has_geometry && (!is.null(axis) || !is.null(coordinate))) {
    stop("`axis`/`coordinate` must not be supplied together with `geometry`.", call. = FALSE)
  }
  if (!has_geometry && (is.null(axis) || is.null(coordinate))) {
    stop("Provide either `geometry`, or both `axis` and `coordinate`.", call. = FALSE)
  }

  if (is.character(image)) image <- ReadImage_fix(image)
  check_sitk_image(image)

  if (!is.null(mask)) {
    if (is.character(mask)) mask <- ReadImage_fix(mask)
    check_sitk_image(mask, arg_name = "mask")
  }

  package_set <- if (has_geometry) {
    .as_slice_package_set(geometry)
  } else {
    build_slice_geometry(image, axis, coordinate)
  }

  list(image = image, mask = mask, package_set = package_set)
}

# Shared per-package/per-slice contour-extraction loop. `binarize_fn`, if
# supplied, is called as binarize_fn(values) -> named list of
# label -> binary (0/1/NaN) vector, and produces one set of contours per
# name at level 0.5 (used by slice_label_contours()); otherwise `levels` is
# contoured directly on the raw sampled values (used by slice_contours()).
.extract_contours <- function(image, mask, package_set, mask_fill,
                               levels = NULL, binarize_fn = NULL) {
  results <- list()

  for (pkg_name in package_set$get_package_names()) {
    pkg <- package_set$get_packages()[[pkg_name]]
    sampled <- pkg$sample_intensity(image, interpolator = "sitkLinear")
    mask_sampled <- if (!is.null(mask)) {
      pkg$sample_intensity(mask, interpolator = "sitkNearestNeighbor")
    } else {
      NULL
    }

    size <- pkg$get_size()
    n_i <- size[1]
    n_j <- size[2]
    spacing <- pkg$get_spacing()
    x_seq <- seq(0, (n_i - 1) * spacing[1], length.out = n_i)
    y_seq <- seq(0, (n_j - 1) * spacing[2], length.out = n_j)

    for (k in sort(unique(sampled$k))) {
      rows <- sampled$k == k
      vals <- sampled$intensity[rows]
      if (!is.null(mask_sampled)) {
        mvals <- mask_sampled$intensity[mask_sampled$k == k]
        vals <- .apply_mask_fill(vals, mvals, mask_fill)
      }

      slice_geom <- pkg$get_slice(k)

      if (is.null(binarize_fn)) {
        z <- .slice_intensity_matrix(vals, n_i, n_j)
        cl <- grDevices::contourLines(x = x_seq, y = y_seq, z = z, levels = levels)
        for (obj_i in seq_along(cl)) {
          path <- cl[[obj_i]]
          world <- .contour_path_to_world(path$x, path$y, slice_geom)
          results[[length(results) + 1]] <- tibble::tibble(
            package = pkg_name, k = k, level = path$level, obj = obj_i,
            vertex = seq_along(path$x),
            x = world[, 1], y = world[, 2], z = world[, 3]
          )
        }
      } else {
        per_label <- binarize_fn(vals)
        for (lv in names(per_label)) {
          z <- .slice_intensity_matrix(per_label[[lv]], n_i, n_j)
          cl <- grDevices::contourLines(x = x_seq, y = y_seq, z = z, levels = 0.5)
          for (obj_i in seq_along(cl)) {
            path <- cl[[obj_i]]
            world <- .contour_path_to_world(path$x, path$y, slice_geom)
            results[[length(results) + 1]] <- tibble::tibble(
              package = pkg_name, k = k, label = as.numeric(lv), obj = obj_i,
              vertex = seq_along(path$x),
              x = world[, 1], y = world[, 2], z = world[, 3]
            )
          }
        }
      }
    }
  }

  results
}

# Bind per-path tibbles from .extract_contours(), drop paths shorter than
# min_vertices (grouped by group_cols), and sort by group_cols + vertex.
# empty_extra_col ("level" or "label") names the empty tibble's extra column
# when results is empty.
.finalize_contours_df <- function(results, group_cols, min_vertices, empty_extra_col) {
  if (length(results) == 0) {
    return(tibble::tibble(
      package = character(), k = integer(),
      !!empty_extra_col := numeric(), obj = integer(), vertex = integer(),
      x = double(), y = double(), z = double()
    ))
  }
  out <- dplyr::bind_rows(results)
  if (!is.null(min_vertices)) {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::filter(dplyr::n() >= min_vertices) %>%
      dplyr::ungroup()
  }
  dplyr::arrange(out, dplyr::across(dplyr::all_of(c(group_cols, "vertex"))))
}

#' Extract iso-intensity contour paths from a slice (or slices) of an image
#'
#' @description
#' The contour-path counterpart to [slice_image()]: instead of per-voxel
#' intensities, returns the vertices of iso-intensity contour lines (via
#' `grDevices::contourLines()`, i.e. marching squares) at one or more
#' requested `levels`, for one or more slices — ready for `ggplot2::geom_path()`
#' rather than a filled raster. Suitable for overlaying a statistical-map
#' threshold or any other iso-intensity boundary on a continuous image.
#'
#' Contouring happens in each slice's own local in-plane coordinates and is
#' then mapped to world coordinates via that slice's `SliceGeometry`, so
#' (unlike the legacy `ggslicer` this replaces) it works for oblique planes,
#' not just axis-aligned ones.
#'
#' @param image A `SimpleITK` image, or a file path (read internally via
#'   [ReadImage_fix()]).
#' @param axis Which image axis is out-of-plane (see [build_slice_geometry()]).
#'   Required unless `geometry` is supplied; must not be supplied together
#'   with `geometry`.
#' @param coordinate One or more world coordinates along `axis` (see
#'   [build_slice_geometry()]). Required unless `geometry` is supplied; must
#'   not be supplied together with `geometry`.
#' @param levels Numeric vector of intensity levels to contour.
#' @param mask An optional `SimpleITK` image or path, sampled
#'   nearest-neighbor onto the same geometry and applied per `mask_fill`
#'   before contouring.
#' @param mask_fill How to treat voxels excluded by `mask` (values `< 0.5`):
#'   `"zero"` (default) zero-fills them, which — since 0 is a value that can
#'   itself be a requested contour level — also traces the mask's own edge
#'   as a contour (useful for drawing a mask/brain outline); `"nan"` fills
#'   them with `NaN` instead, which `contourLines()` skips over without
#'   fabricating a boundary there, better suited to continuous data (e.g. a
#'   signed statistical map) where 0 is itself a meaningful in-range value.
#' @param min_vertices Optional single integer: drop any individual contour
#'   path (a `package`/`k`/`level`/`obj` group) with fewer than this many
#'   vertices, to remove small/noisy paths.
#' @param geometry A pre-built [SliceGeometry]/[SlicePackage]/[SlicePackageSet],
#'   used instead of `axis`/`coordinate`. Must not be supplied together with
#'   `axis`/`coordinate`.
#' @return A tibble with columns `package`, `k`, `level`, `obj` (contour-path
#'   ID within that package/slice/level), `vertex` (1-based order within the
#'   path — always present and pre-sorted, unlike the row-order-dependent
#'   legacy output), `x`, `y`, `z` (world coordinates).
#' @examples
#' \dontrun{
#' image <- ReadImage_fix("statistical_map.nii.gz")
#' mask <- ReadImage_fix("brainmask.nii.gz")
#' contours <- slice_contours(
#'   image, axis = "axial", coordinate = 0, levels = c(-2, 2),
#'   mask = mask, mask_fill = "nan"
#' )
#'
#' library(ggplot2)
#' ggplot(contours, aes(x = x, y = y, group = interaction(level, obj))) +
#'   geom_path()
#' }
#' @export
slice_contours <- function(image, axis = NULL, coordinate = NULL, levels,
                            mask = NULL, mask_fill = c("zero", "nan"),
                            min_vertices = NULL, geometry = NULL) {
  mask_fill <- match.arg(mask_fill)
  if (!is.numeric(levels) || length(levels) < 1) {
    stop("`levels` must be a numeric vector with at least one value.", call. = FALSE)
  }

  inputs <- .resolve_contour_inputs(image, axis, coordinate, mask, geometry)

  results <- .extract_contours(
    inputs$image, inputs$mask, inputs$package_set, mask_fill,
    levels = levels
  )

  .finalize_contours_df(results, c("package", "k", "level", "obj"), min_vertices, "level")
}

#' Extract label-boundary contour paths from a slice (or slices) of a label image
#'
#' @description
#' The discrete-data counterpart to [slice_contours()]: for a label/atlas
#' image, rounds the sampled values to integers, finds the unique nonzero
#' labels present in each slice (optionally restricted to `labels`), and
#' traces each label's own boundary separately (contouring that label's
#' binary indicator at the fixed threshold `0.5`) — ready for
#' `ggplot2::geom_path()` as region-outline annotations, e.g. atlas
#' boundaries drawn over a separately-plotted anatomical image without a
#' solid fill obscuring it.
#'
#' As with [slice_contours()], contouring happens in local in-plane
#' coordinates and is mapped to world coordinates via each slice's own
#' `SliceGeometry`, so it works for oblique planes too.
#'
#' @param image A `SimpleITK` image, or a file path (read internally via
#'   [ReadImage_fix()]), containing integer-valued (or integer-valued-once-
#'   rounded) label data. Always sampled nearest-neighbor, since label data
#'   is discrete.
#' @param axis Which image axis is out-of-plane (see [build_slice_geometry()]).
#'   Required unless `geometry` is supplied; must not be supplied together
#'   with `geometry`.
#' @param coordinate One or more world coordinates along `axis` (see
#'   [build_slice_geometry()]). Required unless `geometry` is supplied; must
#'   not be supplied together with `geometry`.
#' @param labels Optional numeric vector restricting which label values are
#'   contoured. Defaults to every nonzero label value present in each slice.
#' @param mask An optional `SimpleITK` image or path, sampled
#'   nearest-neighbor onto the same geometry and applied per `mask_fill`
#'   before contouring. Since `0` is already the "no label" sentinel that's
#'   always excluded from contouring, `mask_fill` makes no practical
#'   difference here (unlike [slice_contours()]) — it's offered for
#'   consistency with that function's interface, not because the two modes
#'   diverge for label data.
#' @param mask_fill See [slice_contours()]; `"zero"` (default) or `"nan"`.
#' @param min_vertices Optional single integer: drop any individual contour
#'   path (a `package`/`k`/`label`/`obj` group) with fewer than this many
#'   vertices.
#' @param geometry A pre-built [SliceGeometry]/[SlicePackage]/[SlicePackageSet],
#'   used instead of `axis`/`coordinate`. Must not be supplied together with
#'   `axis`/`coordinate`.
#' @return A tibble with columns `package`, `k`, `label`, `obj` (contour-path
#'   ID within that package/slice/label), `vertex` (1-based order within the
#'   path), `x`, `y`, `z` (world coordinates).
#' @examples
#' \dontrun{
#' atlas <- ReadImage_fix("annotation.nii.gz")
#' boundaries <- slice_label_contours(atlas, axis = "axial", coordinate = 0, labels = c(10557, 187))
#'
#' library(ggplot2)
#' ggplot(boundaries, aes(x = x, y = y, group = interaction(label, obj), color = factor(label))) +
#'   geom_path()
#' }
#' @export
slice_label_contours <- function(image, axis = NULL, coordinate = NULL, labels = NULL,
                                  mask = NULL, mask_fill = c("zero", "nan"),
                                  min_vertices = NULL, geometry = NULL) {
  mask_fill <- match.arg(mask_fill)
  if (!is.null(labels) && !is.numeric(labels)) {
    stop("`labels` must be a numeric vector (or NULL).", call. = FALSE)
  }

  inputs <- .resolve_contour_inputs(image, axis, coordinate, mask, geometry)

  binarize_fn <- function(vals) {
    vals <- round(vals)
    present <- sort(unique(vals))
    present <- setdiff(present, 0)
    if (!is.null(labels)) present <- intersect(present, labels)

    out <- list()
    for (lv in present) {
      indicator <- rep(0, length(vals))
      indicator[!is.na(vals) & vals == lv] <- 1
      indicator[is.na(vals)] <- NaN
      out[[as.character(lv)]] <- indicator
    }
    out
  }

  results <- .extract_contours(
    inputs$image, inputs$mask, inputs$package_set, mask_fill,
    binarize_fn = binarize_fn
  )

  .finalize_contours_df(results, c("package", "k", "label", "obj"), min_vertices, "label")
}
