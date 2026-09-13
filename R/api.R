# Coerce a single item to a SlicePackageSet: pass a SlicePackageSet through
# unchanged; wrap a bare SlicePackage or SliceGeometry as a single-package
# set (auto-named "package_1"). Used so sample_images()/slice_image() can
# accept any of the three geometry classes interchangeably as `geometry`.
.as_slice_package_set <- function(x) {
  if (inherits(x, "SlicePackageSet")) {
    return(x)
  }
  if (inherits(x, "SlicePackage") || inherits(x, "SliceGeometry")) {
    return(SlicePackageSet$from_slice_packages(x))
  }
  stop("`geometry` must be a `SliceGeometry`, `SlicePackage`, or `SlicePackageSet` object.", call. = FALSE)
}

#' Names that mark an image as discrete (categorical) data
#'
#' @description
#' The default set of image names that [sample_images()]/[slice_image()]
#' treat as discrete/categorical data (masks, labels, atlases) rather than
#' continuous intensity data, and therefore always sample with
#' `"sitkNearestNeighbor"` regardless of the `interpolator` requested for
#' everything else. Matching is by whole token (splitting the image's name
#' on runs of non-alphanumeric characters), case-insensitively — so
#' `"brain_mask"` matches (token `"mask"`) but `"landmasking_score"` does not
#' (no token equals a discrete name exactly).
#'
#' @return A character vector of names.
#' @export
discrete_data_names <- function() {
  c(
    "mask", "label", "labels", "segmentation", "segmentations", "atlas",
    "seg", "aseg", "aparc", "parcellation", "parcellations", "parcels",
    "roi", "rois", "annotation", "annotations", "regions"
  )
}

# Resolve which interpolator to use for an image given its name: nearest-
# neighbor if any whole token of `name` matches `discrete_names`
# (case-insensitively), otherwise `interpolator` unchanged.
.resolve_interpolator <- function(name, interpolator, discrete_names) {
  tokens <- strsplit(tolower(name), "[^A-Za-z0-9]+")[[1]]
  tokens <- tokens[nzchar(tokens)]
  if (any(tokens %in% tolower(discrete_names))) {
    return("sitkNearestNeighbor")
  }
  interpolator
}

#' Sample one or more SimpleITK images onto a shared slice geometry
#'
#' @description
#' The core sampling primitive behind [build_slice_geometry()]/[slice_image()]:
#' given a geometry (any of [SliceGeometry], [SlicePackage], or
#' [SlicePackageSet] — auto-coerced to a `SlicePackageSet`) and one or more
#' named images, returns a single tidy data frame with one row per sample
#' point and one column per image.
#'
#' Each image may be given as an already-loaded `SimpleITK` image or as a
#' file path (read internally via [ReadImage_fix()]). The interpolator used
#' for each image is `interpolator` by default, except that images whose
#' name matches [discrete_data_names()] (e.g. `"mask"`, `"label"`, `"atlas"`)
#' always use `"sitkNearestNeighbor"`, and `interpolator_overrides` always
#' wins over both. Images with more than 3 dimensions are indexed via
#' `extra_index`; images that are exactly 3D ignore `extra_index` and have
#' their (single) sampled value broadcast across every combination present
#' from higher-dimensional images.
#'
#' @param geometry A [SliceGeometry], [SlicePackage], or [SlicePackageSet] object.
#' @param images A single `SimpleITK` image or file path (sampled into an
#'   `"intensity"` column), or a fully named list of images/paths (sampled
#'   into columns named after the list).
#' @param interpolator Default SimpleITK interpolator name (e.g.
#'   `"sitkLinear"`, the default, or `"sitkNearestNeighbor"`), used for any
#'   image not matched by `discrete_names` or `interpolator_overrides`.
#' @param interpolator_overrides A named list mapping an image's name to an
#'   explicit interpolator, taking precedence over both `interpolator` and
#'   the `discrete_names` rule.
#' @param discrete_names Character vector of names (matched by whole token,
#'   case-insensitively) that always use `"sitkNearestNeighbor"`. Defaults to
#'   [discrete_data_names()].
#' @param extra_index A named list of 0-indexed index vectors, one per
#'   non-spatial dimension shared by any higher-dimensional (4D/5D) image in
#'   `images` (e.g. `list(t = 0:9)`). Images that are exactly 3D ignore this
#'   and are sampled once, broadcasting across every combination.
#' @return A tibble with columns `package`, any names in `extra_index`, `i`,
#'   `j`, `k`, `x`, `y`, `z`, and one column per name in `images`.
#' @export
sample_images <- function(geometry, images, interpolator = "sitkLinear",
                           interpolator_overrides = list(),
                           discrete_names = discrete_data_names(),
                           extra_index = list()) {
  pkgset <- .as_slice_package_set(geometry)

  if (!is.list(images)) {
    images <- list(intensity = images)
  } else if (is.null(names(images)) || any(!nzchar(names(images)))) {
    stop("`images` must be a fully named list when supplying more than one image.", call. = FALSE)
  }

  n_extra <- length(extra_index)
  base <- pkgset$get_sample_points()
  if (n_extra == 0) {
    out <- base
  } else {
    combos <- expand.grid(extra_index, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    out <- dplyr::bind_rows(lapply(seq_len(nrow(combos)), function(r) {
      rows <- base
      for (nm in names(extra_index)) rows[[nm]] <- combos[[nm]][r]
      rows
    }))
    out <- dplyr::select(out, "package", dplyr::all_of(names(extra_index)), dplyr::everything())
  }

  for (name in names(images)) {
    img <- images[[name]]
    if (is.character(img)) img <- ReadImage_fix(img)
    check_sitk_image(img, arg_name = paste0("images$", name))

    dim <- img$GetDimension()
    needed_n_extra <- dim - 3
    if (needed_n_extra == 0) {
      img_extra_index <- list()
    } else if (needed_n_extra == n_extra) {
      img_extra_index <- extra_index
    } else {
      stop(sprintf(
        "Image `%s` has dimension %d, which needs `extra_index` with %d entries, but `extra_index` has %d.",
        name, dim, needed_n_extra, n_extra
      ), call. = FALSE)
    }

    interp <- interpolator_overrides[[name]]
    if (is.null(interp)) {
      interp <- .resolve_interpolator(name, interpolator, discrete_names)
    }

    sampled <- pkgset$sample_intensity(img, extra_index = img_extra_index, interpolator = interp)
    join_cols <- c("package", names(img_extra_index), "i", "j", "k")
    sampled <- sampled[, c(join_cols, "intensity")]
    colnames(sampled)[colnames(sampled) == "intensity"] <- name

    out <- dplyr::left_join(out, sampled, by = join_cols)
  }

  out
}

#' Build a slice geometry from an image and world coordinates along one axis
#'
#' @description
#' The "smart" geometry builder behind [slice_image()]'s default (non-escape-
#' hatch) path: given an image, an out-of-plane axis, and one or more world
#' coordinates along that axis, builds the corresponding axis-aligned
#' slice(s) (via [SliceGeometry$from_image_axis()][SliceGeometry]) and always
#' returns a [SlicePackageSet], for a uniform return type regardless of input:
#'
#' - A single `coordinate` gives a `SlicePackageSet` with one single-slice package.
#' - Multiple, evenly-spaced coordinates give a `SlicePackageSet` with one
#'   regular [SlicePackage] (so the whole stack resamples in a single call).
#' - Multiple, unevenly-spaced coordinates give a `SlicePackageSet` with one
#'   single-slice package per coordinate (since a [SlicePackage] cannot
#'   represent uneven spacing) — auto-detected, not something the caller
#'   needs to specify.
#'
#' @param image A 3D `SimpleITK` image.
#' @param axis Which image axis is out-of-plane: `1`/`2`/`3`, `"x"`/`"y"`/`"z"`,
#'   or (assuming right-anterior-superior orientation)
#'   `"sagittal"`/`"coronal"`/`"axial"`/`"horizontal"`.
#' @param coordinate A numeric vector of one or more world coordinates along
#'   `axis` (each snapped to the nearest voxel plane, as in
#'   `SliceGeometry$from_image_axis()`). Need not be evenly spaced.
#' @return A new [SlicePackageSet].
#' @export
build_slice_geometry <- function(image, axis, coordinate) {
  check_sitk_image(image)
  if (!is.numeric(coordinate) || length(coordinate) < 1) {
    stop("`coordinate` must be a numeric vector with at least one value.", call. = FALSE)
  }

  slices <- lapply(coordinate, function(co) SliceGeometry$from_image_axis(image, axis, co))

  if (length(slices) == 1) {
    return(SlicePackageSet$from_slice_packages(slices[[1]]))
  }

  pkg <- tryCatch(SlicePackage$from_slices(slices), error = function(e) NULL)
  if (!is.null(pkg)) {
    return(SlicePackageSet$from_slice_packages(pkg))
  }
  SlicePackageSet$from_slice_packages(slices)
}

#' Sample a slice (or slices) of an image, ready for grammar-of-graphics plotting
#'
#' @description
#' The main user-facing entry point: given an image, an axis, and one or more
#' world coordinates along it (or, via the `geometry` escape hatch, a
#' hand-built [SliceGeometry]/[SlicePackage]/[SlicePackageSet] for oblique or
#' custom cases), returns a tidy data frame with one row per sample point and
#' one column per image (the main `image`, named `"value"`, plus any
#' `extra_images`) — ready to hand to `ggplot2`.
#'
#' Internally, this is [build_slice_geometry()] (unless `geometry` is
#' supplied directly) followed by [sample_images()]; see those for the
#' geometry-construction and interpolator-selection details.
#'
#' @param image A `SimpleITK` image, or a file path (read internally via
#'   [ReadImage_fix()]).
#' @param axis Which image axis is out-of-plane (see [build_slice_geometry()]).
#'   Required unless `geometry` is supplied; must not be supplied together
#'   with `geometry`.
#' @param coordinate One or more world coordinates along `axis` (see
#'   [build_slice_geometry()]). Required unless `geometry` is supplied; must
#'   not be supplied together with `geometry`.
#' @param extra_images A fully named list of additional images/paths to
#'   sample onto the same geometry (e.g. `list(mask = "brainmask.nii.gz")`),
#'   sampled into columns named after the list.
#' @param interpolator Default SimpleITK interpolator name for `image` and
#'   any `extra_images` not otherwise matched; see [sample_images()].
#' @param interpolator_overrides A named list of explicit per-image
#'   interpolator overrides; see [sample_images()].
#' @param discrete_names Character vector of names treated as discrete
#'   (categorical) data; see [discrete_data_names()].
#' @param extra_index Non-spatial index selection for 4D/5D images; see
#'   [sample_images()].
#' @param geometry A pre-built [SliceGeometry]/[SlicePackage]/[SlicePackageSet],
#'   used instead of `axis`/`coordinate` (e.g. for an oblique DICOM-derived
#'   geometry). Must not be supplied together with `axis`/`coordinate`.
#' @return A tibble; see [sample_images()].
#' @export
slice_image <- function(image, axis = NULL, coordinate = NULL, extra_images = list(),
                         interpolator = "sitkLinear", interpolator_overrides = list(),
                         discrete_names = discrete_data_names(), extra_index = list(),
                         geometry = NULL) {
  has_geometry <- !is.null(geometry)
  if (has_geometry && (!is.null(axis) || !is.null(coordinate))) {
    stop("`axis`/`coordinate` must not be supplied together with `geometry`.", call. = FALSE)
  }
  if (!has_geometry && (is.null(axis) || is.null(coordinate))) {
    stop("Provide either `geometry`, or both `axis` and `coordinate`.", call. = FALSE)
  }

  if (is.character(image)) image <- ReadImage_fix(image)
  check_sitk_image(image)

  geom <- if (has_geometry) geometry else build_slice_geometry(image, axis, coordinate)

  images <- c(list(value = image), extra_images)

  sample_images(
    geom, images,
    interpolator = interpolator,
    interpolator_overrides = interpolator_overrides,
    discrete_names = discrete_names,
    extra_index = extra_index
  )
}
