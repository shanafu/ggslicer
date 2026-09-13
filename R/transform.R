# Read a file path as a SimpleITK transform, routing MINC .xfm files through
# read_minc_transform() (SimpleITK::ReadTransform() silently mis-parses any
# .xfm containing a Grid_Transform block or multiple concatenated blocks --
# confirmed directly, not assumed; see CLAUDE.md) and everything else through
# the normal SimpleITK reader.
.read_transform_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "xfm") {
    read_minc_transform(path)
  } else {
    SimpleITK::ReadTransform(path)
  }
}

#' Read a MINC transform (`.xfm`) file as a SimpleITK transform
#'
#' @description
#' Parses an MNI transform file directly (no dependency beyond what this
#' package already requires): one or more `Transform_Type = Linear;` blocks
#' (a 3x4 matrix) and/or `Transform_Type = Grid_Transform;` blocks (a
#' reference to a companion displacement-field MINC volume, resolved
#' relative to `path`'s own directory), optionally concatenated in one file.
#' An `Invert_Flag` on a `Grid_Transform` block is ignored — confirmed
#' (round-tripping a real forward/inverse transform pair) that the
#' referenced displacement volume already contains the correct-direction
#' field, needing no extra sign handling.
#'
#' `SimpleITK::ReadTransform()` must not be used for a `.xfm` file containing
#' a `Grid_Transform` block or multiple concatenated blocks — confirmed it
#' silently returns a transform with the wrong dimensions and all-zero
#' displacement rather than erroring. (A pure single-`Linear`-block `.xfm`
#' does parse correctly via `SimpleITK::ReadTransform()`, but this function
#' handles that case too, for a single entry point regardless of content.)
#'
#' @param path Path to a `.xfm` file.
#' @return A `SimpleITK` transform: a single `AffineTransform` or
#'   `DisplacementFieldTransform` if `path` has exactly one block, or a
#'   `CompositeTransform` (applying the blocks in file order) if it has more
#'   than one.
#' @examples
#' \dontrun{
#' xfm <- read_minc_transform("registration.xfm")
#' xfm$TransformPoint(c(0, 0, 0))
#' }
#' @export
read_minc_transform <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!grepl("MNI Transform File", text, fixed = TRUE)) {
    stop("`path` does not look like an MNI transform file (.xfm): ", path, call. = FALSE)
  }

  type_matches <- gregexpr("Transform_Type\\s*=\\s*(\\w+)\\s*;", text, perl = TRUE)[[1]]
  if (type_matches[1] == -1) {
    stop("No `Transform_Type` blocks found in: ", path, call. = FALSE)
  }
  match_lengths <- attr(type_matches, "match.length")
  type_names <- vapply(seq_along(type_matches), function(i) {
    block_header <- substr(text, type_matches[i], type_matches[i] + match_lengths[i] - 1)
    sub("Transform_Type\\s*=\\s*(\\w+)\\s*;", "\\1", block_header, perl = TRUE)
  }, character(1))

  block_starts <- type_matches + match_lengths
  block_ends <- c(type_matches[-1] - 1, nchar(text))
  base_dir <- dirname(path)

  transforms <- vector("list", length(type_names))
  for (i in seq_along(type_names)) {
    ttype <- type_names[i]
    block_text <- substr(text, block_starts[i], block_ends[i])

    if (ttype == "Linear") {
      m <- regmatches(block_text, regexpr("Linear_Transform\\s*=([^;]+);", block_text, perl = TRUE))
      if (length(m) == 0 || !nzchar(m)) {
        stop("Malformed `Linear_Transform` block in: ", path, call. = FALSE)
      }
      nums_text <- sub("^Linear_Transform\\s*=", "", m)
      nums_text <- sub(";\\s*$", "", nums_text)
      nums <- suppressWarnings(as.numeric(strsplit(trimws(nums_text), "\\s+")[[1]]))
      if (length(nums) != 12 || anyNA(nums)) {
        stop("Malformed `Linear_Transform` block (expected 12 numbers) in: ", path, call. = FALSE)
      }
      mat <- matrix(nums, nrow = 3, byrow = TRUE)
      t <- SimpleITK::AffineTransform(3L)
      t$SetMatrix(as.vector(t(mat[, 1:3])))
      t$SetTranslation(mat[, 4])
      transforms[[i]] <- t
    } else if (ttype == "Grid_Transform") {
      m <- regmatches(block_text, regexpr("Displacement_Volume\\s*=\\s*([^;]+);", block_text, perl = TRUE))
      if (length(m) == 0 || !nzchar(m)) {
        stop("Malformed `Grid_Transform` block in: ", path, call. = FALSE)
      }
      vol_path <- sub("^Displacement_Volume\\s*=\\s*", "", m)
      vol_path <- trimws(sub(";\\s*$", "", vol_path))
      if (!grepl("^(/|[A-Za-z]:)", vol_path)) {
        vol_path <- file.path(base_dir, vol_path)
      }
      # Read plainly -- never through ReadImage_fix()/orientation_correction().
      grid_img <- SimpleITK::ReadImage(vol_path, "sitkVectorFloat64")
      transforms[[i]] <- SimpleITK::DisplacementFieldTransform(grid_img)
    } else {
      stop(
        "Unsupported `Transform_Type` \"", ttype, "\" in ", path,
        " (only \"Linear\" and \"Grid_Transform\" are supported).",
        call. = FALSE
      )
    }
  }

  if (length(transforms) == 1) {
    return(transforms[[1]])
  }

  # File blocks are meant to apply in file order [T1, T2, ...]; SimpleITK's
  # CompositeTransform applies the *last*-added transform first, so add them
  # in reverse (verified with a non-commuting synthetic case; see CLAUDE.md).
  composite <- SimpleITK::CompositeTransform(3L)
  for (t in rev(transforms)) {
    composite$AddTransform(t)
  }
  composite
}

#' Apply a SimpleITK transform to a tidy data frame's world coordinates
#'
#' @description
#' A small, generic utility: given any tidy data frame with `x`/`y`/`z`
#' (world coordinate) columns — [slice_grid()]'s output, [slice_contours()]'s
#' output, or anything else — replaces those columns with their positions
#' after applying `transform`, leaving every other column untouched. This is
#' what makes visualizing a registration warp possible: build a grid (or
#' contour, or any point set) on the fixed/reference image, then transform
#' its points through the registration transform.
#'
#' Transforms one point at a time (`SimpleITK` has no vectorized/batch
#' transform-point API), so this can take a while for a large data frame —
#' progress is reported via `pbapply::pbapply()`.
#'
#' @param df A data frame with `x_col`/`y_col`/`z_col` columns.
#' @param transform A `SimpleITK` transform object, or a file path (read via
#'   [read_minc_transform()] for `.xfm`, or `SimpleITK::ReadTransform()`
#'   otherwise).
#' @param invert If `TRUE`, apply the inverse of `transform`
#'   (`transform$GetInverse()`). Works for affine/linear transforms; SimpleITK
#'   does not support inverting a `DisplacementFieldTransform` (or a
#'   composite containing one) this way and will error clearly if asked —
#'   load the separately-computed inverse-warp file instead (the standard
#'   registration-tool convention, and why e.g. ANTs always writes both
#'   `*Warp.nii.gz` and `*InverseWarp.nii.gz`).
#' @param x_col,y_col,z_col Names of the world-coordinate columns to
#'   transform.
#' @return `df`, with `x_col`/`y_col`/`z_col` replaced by their transformed
#'   coordinates.
#' @examples
#' \dontrun{
#' image <- ReadImage_fix("fixed_template.nii")
#' grid_df <- slice_grid(image, axis = "axial", coordinate = 0, spacing = 5)
#'
#' # Visualize a registration warp by transforming a grid built on the fixed image
#' warped_df <- transform_points(grid_df, "registration.xfm")
#' }
#' @export
transform_points <- function(df, transform, invert = FALSE, x_col = "x", y_col = "y", z_col = "z") {
  if (is.character(transform)) {
    transform <- .read_transform_file(transform)
  }
  if (invert) {
    transform <- transform$GetInverse()
  }

  coords <- cbind(df[[x_col]], df[[y_col]], df[[z_col]])
  transformed <- pbapply::pbapply(coords, 1, function(point) transform$TransformPoint(point))

  df[[x_col]] <- transformed[1, ]
  df[[y_col]] <- transformed[2, ]
  df[[z_col]] <- transformed[3, ]
  df
}
