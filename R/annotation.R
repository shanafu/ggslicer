# The area-weighted centroid of a closed, planar polygon embedded in 3D
# (vertices already form a closed loop -- confirmed directly that both
# grDevices::contourLines() and contourpy always return `x[1] == x[n]`,
# `y[1] == y[n]` -- so no separate wraparound edge is needed here; consecutive
# vertex pairs already close the loop). Computed via a triangle fan from the
# vertex mean (used purely for numerical conditioning, not as a geometric
# assumption -- the result is exact regardless of that reference point's
# position, including for concave shapes) using 3D cross products, so it
# needs no local 2D (i, j) reparameterization -- important because
# contour_centroids() only ever sees world x/y/z coordinates, with no
# SliceGeometry available to recover a local in-plane basis from.
# Returns list(x=, y=, z=, area=). Falls back to the plain vertex mean when
# the enclosed area is ~0 (a degenerate/collinear contour).
.polygon_centroid_3d <- function(x, y, z) {
  v <- cbind(x, y, z)
  n <- nrow(v)
  ref <- colMeans(v)
  a <- sweep(v, 2, ref, "-")

  edge_a <- a[-n, , drop = FALSE]
  edge_b <- a[-1, , drop = FALSE]

  cross_mat <- cbind(
    edge_a[, 2] * edge_b[, 3] - edge_a[, 3] * edge_b[, 2],
    edge_a[, 3] * edge_b[, 1] - edge_a[, 1] * edge_b[, 3],
    edge_a[, 1] * edge_b[, 2] - edge_a[, 2] * edge_b[, 1]
  )
  tri_area_vec <- 0.5 * cross_mat
  tri_centroid_rel <- (edge_a + edge_b) / 3

  total_area_vec <- colSums(tri_area_vec)
  total_area <- sqrt(sum(total_area_vec^2))

  if (total_area < sqrt(.Machine$double.eps)) {
    return(list(x = ref[1], y = ref[2], z = ref[3], area = 0))
  }

  n_hat <- total_area_vec / total_area
  s <- as.numeric(tri_area_vec %*% n_hat)
  centroid_rel <- colSums(s * tri_centroid_rel) / sum(s)
  centroid <- ref + centroid_rel

  list(x = centroid[1], y = centroid[2], z = centroid[3], area = total_area)
}

# Look up a display name for `label` from `label_names` (a named vector,
# names matched as the integer-rounded label value), falling back to the
# label's own integer value as a string when unmapped or when
# `label_names` is NULL. No lookup table for label -> region name ships
# with this package or with `testdata/` (confirmed directly -- real atlases
# like DSURQE keep this in a separate CSV, not in the label image itself),
# so this must come from the caller.
.label_display_name <- function(label, label_names) {
  key <- as.character(as.integer(round(label)))
  if (!is.null(label_names) && key %in% names(label_names)) {
    unname(label_names[[key]])
  } else {
    key
  }
}

#' Compute one label per connected component from a slice_label_contours() data frame
#'
#' @description
#' Given the output of [slice_label_contours()], finds the centroid of each
#' connected component (each `package`/`k`/`label`/`obj` group already
#' distinguishes separate components sharing the same label -- a label
#' present as multiple disjoint regions in one slice gets one row per
#' region here, not one row per label) -- ready for
#' `ggplot2::geom_text(aes(x = x, y = y, label = name))` as region-name
#' annotations.
#'
#' The centroid is the true, area-weighted centroid of the enclosed polygon
#' (a triangle-fan computation, done directly in world coordinates), not
#' just the mean of the traced boundary's vertices -- more robust for
#' irregular/concave region shapes,
#' though even this can occasionally fall outside a very non-convex (e.g.
#' crescent-shaped) region; inspect placements before relying on them for a
#' specific figure.
#'
#' @param contours_df A data frame as returned by [slice_label_contours()]
#'   (columns `package`, `k`, `label`, `obj`, `vertex`, `x`, `y`, `z`).
#' @param label_names Optional named vector mapping a label's integer value
#'   (as a string name, e.g. `c("187" = "Hippocampus")`) to its display
#'   text. Labels not present in `label_names` (or when `label_names` is
#'   `NULL`) fall back to their own integer value as a string.
#' @param min_area Optional single number: drop any connected component
#'   whose true polygon area (world-space units², e.g. mm^2) is below this
#'   value -- a more accurate size filter than
#'   [slice_label_contours()]'s own `min_vertices` (which really tracks
#'   perimeter/boundary complexity, not enclosed area).
#' @return A tibble with columns `package`, `k`, `label`, `name` (display
#'   text), `obj`, `x`, `y`, `z` (centroid world coordinates), `area`
#'   (world-space units²).
#' @examples
#' \dontrun{
#' atlas <- ReadImage_fix("atlas_labels.nii.gz")
#' contours <- slice_label_contours(atlas, axis = "axial", coordinate = 0)
#' annotations <- contour_centroids(
#'   contours,
#'   label_names = c("187" = "Hippocampus", "10557" = "Cortex"),
#'   min_area = 2
#' )
#'
#' library(ggplot2)
#' ggplot(annotations, aes(x = x, y = y, label = name)) + geom_text()
#' }
#' @export
contour_centroids <- function(contours_df, label_names = NULL, min_area = NULL) {
  required_cols <- c("package", "k", "label", "obj", "vertex", "x", "y", "z")
  if (!all(required_cols %in% names(contours_df))) {
    stop(
      "`contours_df` must have columns ", paste(required_cols, collapse = ", "),
      " (as returned by slice_label_contours()).",
      call. = FALSE
    )
  }
  if (!is.null(min_area) && (!is.numeric(min_area) || length(min_area) != 1 || min_area < 0)) {
    stop("`min_area` must be a single non-negative number (or NULL).", call. = FALSE)
  }

  if (nrow(contours_df) == 0) {
    return(tibble::tibble(
      package = character(), k = integer(), label = numeric(), name = character(),
      obj = integer(), x = double(), y = double(), z = double(), area = double()
    ))
  }

  ordered <- dplyr::arrange(contours_df, .data$package, .data$k, .data$label, .data$obj, .data$vertex)
  groups <- dplyr::group_by(ordered, .data$package, .data$k, .data$label, .data$obj)

  out <- dplyr::summarise(
    groups,
    .centroid = list(.polygon_centroid_3d(.data$x, .data$y, .data$z)),
    .groups = "drop"
  )
  out$x <- vapply(out$.centroid, function(c) c$x, numeric(1))
  out$y <- vapply(out$.centroid, function(c) c$y, numeric(1))
  out$z <- vapply(out$.centroid, function(c) c$z, numeric(1))
  out$area <- vapply(out$.centroid, function(c) c$area, numeric(1))
  out$.centroid <- NULL

  if (!is.null(min_area)) {
    out <- out[out$area >= min_area, ]
  }

  out$name <- vapply(out$label, .label_display_name, character(1), label_names = label_names)
  dplyr::select(out, "package", "k", "label", "name", "obj", "x", "y", "z", "area")
}

#' Place region-name text at the centroid of each labeled region in a slice
#'
#' @description
#' A one-call convenience wrapper: runs [slice_label_contours()] then
#' [contour_centroids()] on the result -- for placing region-name text
#' annotations directly from an image/axis/coordinate, without first
#' building the contours data frame yourself. If you already need
#' `slice_label_contours()`'s output separately (e.g. to also draw region
#' outlines via `geom_path()`), call [contour_centroids()] directly on that
#' existing result instead, to avoid recomputing the same contours twice.
#'
#' @inheritParams slice_label_contours
#' @param label_names Optional named vector mapping a label's integer value
#'   (as a string name) to its display text; see [contour_centroids()].
#' @param min_area Optional single number: drop any connected component
#'   whose true polygon area is below this value; see [contour_centroids()].
#'   Independent of `min_vertices`, which is [slice_label_contours()]'s own
#'   perimeter/vertex-count-based filter.
#' @return A tibble; see [contour_centroids()].
#' @examples
#' \dontrun{
#' atlas <- ReadImage_fix("atlas_labels.nii.gz")
#' annotations <- slice_label_annotations(
#'   atlas, axis = "axial", coordinate = 0,
#'   label_names = c("187" = "Hippocampus", "10557" = "Cortex"),
#'   min_area = 2
#' )
#'
#' library(ggplot2)
#' ggplot(annotations, aes(x = x, y = y, label = name)) + geom_text()
#' }
#' @export
slice_label_annotations <- function(image, axis = NULL, coordinate = NULL, labels = NULL,
                                     label_names = NULL,
                                     mask = NULL, mask_fill = c("zero", "nan"),
                                     min_vertices = NULL, min_area = NULL, geometry = NULL) {
  mask_fill <- match.arg(mask_fill)
  contours_df <- slice_label_contours(
    image, axis = axis, coordinate = coordinate, labels = labels,
    mask = mask, mask_fill = mask_fill, min_vertices = min_vertices, geometry = geometry
  )
  contour_centroids(contours_df, label_names = label_names, min_area = min_area)
}

#' Ready-made ggplot2 layer for region-name annotations produced by slice_label_annotations()
#'
#' @description
#' Wraps `ggplot2::geom_text()` for a [slice_label_annotations()]/
#' [contour_centroids()] data frame, so it can be added to a plot directly --
#' the same ready-made-layer convenience `slice_grid_layers()`/
#' `slice_warp_arrows_layer()` provide for their own outputs.
#'
#' @param annotations_df A data frame as returned by
#'   [slice_label_annotations()]/[contour_centroids()].
#' @param color,size,alpha,fontface Style for the text labels, passed to
#'   `ggplot2::geom_text()`.
#' @return A single `ggplot2::geom_text()` layer. Add it (e.g. `plt + layer`)
#'   to a plot that already maps `x`/`y` via its own `aes()`.
#' @examples
#' \dontrun{
#' atlas <- ReadImage_fix("atlas_labels.nii.gz")
#' annotations <- slice_label_annotations(atlas, axis = "axial", coordinate = 0)
#' layer <- slice_label_annotations_layer(annotations)
#'
#' library(ggplot2)
#' ggplot(annotations, aes(x = x, y = y)) + layer
#' }
#' @importFrom rlang .data
#' @export
slice_label_annotations_layer <- function(annotations_df, color = "black", size = 3.5,
                                           alpha = 1, fontface = "plain") {
  ggplot2::geom_text(
    data = annotations_df,
    mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
    color = color, size = size, alpha = alpha, fontface = fontface
  )
}
