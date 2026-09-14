#' Validate that `image` is a SimpleITK image object
#'
#' Gives a specifically helpful error if the caller passed a file path
#' instead (an easy mistake: several functions elsewhere in this package,
#' e.g. `ReadImage_fix()`, take a path and read it internally, but the
#' functions that call this validator expect an already-loaded image).
#'
#' @param image The value to check.
#' @param arg_name Name to use for `image` in the error message.
#'
#' @return `TRUE`, invisibly, if `image` is a valid SimpleITK image; otherwise errors.
#' @examples
#' \dontrun{
#' image <- ReadImage_fix("brain_image.nii")
#' check_sitk_image(image) # TRUE, invisibly
#' check_sitk_image("brain_image.nii") # errors with a clear message
#' }
#' @keywords internal
check_sitk_image <- function(image, arg_name = "image") {
  if (inherits(image, "_p_itk__simple__Image")) {
    return(invisible(TRUE))
  }
  if (is.character(image)) {
    stop(
      "`", arg_name, "` must be a SimpleITK image object, not a file path ",
      "(got \"", image[1], "\"). Read the file first, e.g. `",
      arg_name, " <- ReadImage_fix(\"", image[1], "\")`, then pass that.",
      call. = FALSE
    )
  }
  stop(
    "`", arg_name, "` must be a SimpleITK image object (e.g. from `ReadImage_fix()` ",
    "or `SimpleITK::ReadImage()`), not a ", paste(class(image), collapse = "/"), ".",
    call. = FALSE
  )
}

# Whether any whole token of `name` (split on runs of non-alphanumeric
# characters, case-insensitively) matches one of `discrete_names` exactly.
# Shared by api.R's .resolve_interpolator() and suggest_contour_levels(),
# so the token-matching rule only lives in one place.
.matches_discrete_name <- function(name, discrete_names) {
  tokens <- strsplit(tolower(name), "[^A-Za-z0-9]+")[[1]]
  tokens <- tokens[nzchar(tokens)]
  any(tokens %in% tolower(discrete_names))
}

# Whether `values` looks like discrete/categorical data: mostly integer-
# valued and relatively few distinct values. A data-driven safety net for
# inputs with no name to check at all (a raw image) or an unusually-named
# discrete column -- used alongside, not instead of, .matches_discrete_name().
# Thresholds are a reasonable starting heuristic, not empirically tuned.
.looks_discrete <- function(values, integer_frac_threshold = 0.99, max_unique = 50) {
  v <- values[!is.na(values)]
  if (length(v) == 0) {
    return(FALSE)
  }
  frac_integer <- mean(abs(v - round(v)) < 1e-6)
  n_unique <- length(unique(v))
  frac_integer >= integer_frac_threshold && n_unique <= max_unique
}

# n interior quantiles of values (excludes the 0th/100th percentile, since
# the min/max aren't useful contour levels).
.quantile_levels <- function(values, n) {
  values <- values[!is.na(values)]
  probs <- seq(0, 1, length.out = n + 2)[-c(1, n + 2)]
  unname(stats::quantile(values, probs = probs))
}

# At least min_n local minima ("troughs") of a kernel density estimate of
# values -- natural boundaries between distinct populations/tissue classes.
# Retries with a progressively narrower bandwidth if too few troughs are
# found (a unimodal distribution has none at any bandwidth), then pads out
# any remaining shortfall with quantile levels.
.trough_levels <- function(values, min_n) {
  values <- values[!is.na(values)]
  range_width <- diff(range(values))
  if (range_width <= 0) {
    return(rep(values[1], min_n))
  }

  # A local minimum only counts as a real "trough" if it's a genuine dip
  # relative to a visibly higher point on *both* sides (a topographic-
  # prominence-style filter), not just numerical wiggle in the density
  # estimate -- e.g. a truly unimodal distribution's KDE can still show tiny
  # spurious local minima out in its low-density tails, which have ~zero
  # prominence and must not be reported as troughs.
  find_troughs <- function(bw, prominence_frac = 0.02) {
    d <- stats::density(values, bw = bw)
    y <- d$y
    dy <- diff(y)
    sgn <- sign(dy)
    local_min_idx <- which(diff(sgn) == 2) + 1
    if (length(local_min_idx) == 0) {
      return(numeric(0))
    }
    local_max_idx <- which(diff(sgn) == -2) + 1
    boundary_idx <- c(1, length(y))
    extrema_idx <- sort(c(local_max_idx, boundary_idx))
    threshold <- prominence_frac * max(y)

    keep <- vapply(local_min_idx, function(idx) {
      left <- extrema_idx[extrema_idx < idx]
      right <- extrema_idx[extrema_idx > idx]
      left_max <- if (length(left) > 0) max(y[max(left)]) else y[idx]
      right_max <- if (length(right) > 0) max(y[min(right)]) else y[idx]
      (min(left_max, right_max) - y[idx]) >= threshold
    }, logical(1))

    d$x[local_min_idx[keep]]
  }

  # Deliberately *oversmoothed*, and deliberately *fixed* (not retried at a
  # narrower bandwidth if too few troughs are found): bw.nrd0() alone (a
  # plug-in rule of thumb tuned for overall density-estimate accuracy, not a
  # clean mode count) is already narrow enough to produce spurious noise
  # troughs alongside/instead of the real one (bw.SJ() is worse still --
  # confirmed ~6x narrower than bw.nrd0() on a test bimodal sample, which
  # reproduced exactly this problem). Retrying at a narrower bandwidth
  # whenever min_n isn't met was tried and confirmed empirically *not* to
  # help -- for a genuinely unimodal distribution (zero real troughs at any
  # bandwidth), narrowing only ever surfaces sampling-noise "troughs" that
  # still pass the prominence filter, never a real one; simply padding the
  # shortfall with quantile levels (below) is safer and no less accurate.
  bw <- 2 * stats::bw.nrd0(values)
  troughs <- find_troughs(bw)

  if (length(troughs) < min_n) {
    extra <- min_n - length(troughs)
    troughs <- c(troughs, .quantile_levels(values, extra))
  }
  sort(troughs)
}

#' Suggest a set of contour levels from an image or a sampled data frame
#'
#' @description
#' Computes a reasonable, automatic set of intensity levels to pass as
#' [slice_contours()]'s `levels` argument, instead of picking them by hand.
#' Works on either a whole image (using its full voxel-intensity
#' distribution) or a [slice_image()]/[sample_images()]-output data frame
#' (using one named column) -- both reduce to one plain numeric vector, and
#' everything past that point is shared.
#'
#' Two methods are available. `"quantile"` (the default) returns `n`
#' evenly-spaced interior quantiles — simple, deterministic, and always
#' returns exactly `n` levels regardless of the data's shape. `"troughs"`
#' instead finds local minima of a kernel density estimate — the natural
#' boundaries between distinct populations in the data (e.g. tissue classes
#' in an anatomical image) — returning at least `min_n` of them; if fewer
#' troughs exist than `min_n` (a unimodal distribution has none at all), the
#' shortfall is padded out with quantile levels. Because this method depends
#' on bandwidth selection, results are a reasonable heuristic, not a
#' guaranteed-optimal set of levels — inspect them before trusting them for
#' a specific analysis claim.
#'
#' `levels` computed this way don't make sense for discrete/categorical data
#' (masks, labels, atlases) — this errors clearly if it looks like `x`
#' (or the selected `column`) is one, checked two ways: by name (matching
#' [discrete_data_names()], exactly like [sample_images()]'s interpolator
#' selection) and, since a raw image has no name to check, by inspecting the
#' values themselves (mostly integer-valued, relatively few distinct values).
#' Use [slice_label_contours()] for that data instead.
#'
#' @param x A `SimpleITK` image, a file path (read internally via
#'   [ReadImage_fix()]), or a data frame (e.g. from [slice_image()]).
#' @param column When `x` is a data frame, the column to use. Ignored
#'   otherwise. Default `"value"` matches [slice_image()]'s main-image column.
#' @param method `"quantile"` (default) or `"troughs"`.
#' @param n Number of levels to return, for `method = "quantile"`.
#' @param min_n Minimum number of levels to return, for `method = "troughs"`.
#' @param discrete_names Character vector of names (matched by whole token,
#'   case-insensitively) that mark `column` as discrete data. Defaults to
#'   [discrete_data_names()].
#' @return A numeric vector of suggested contour levels.
#' @examples
#' \dontrun{
#' image <- ReadImage_fix("brain_image.nii")
#'
#' # Directly from the whole image
#' levels <- suggest_contour_levels(image, method = "quantile", n = 5)
#'
#' # From a slice_image() data frame instead
#' df <- slice_image(image, axis = "axial", coordinate = 0)
#' levels <- suggest_contour_levels(df, column = "value", method = "troughs", min_n = 2)
#'
#' contours <- slice_contours(image, axis = "axial", coordinate = 0, levels = levels)
#' }
#' @export
suggest_contour_levels <- function(x, column = "value", method = c("quantile", "troughs"),
                                    n = 5, min_n = 1, discrete_names = discrete_data_names()) {
  method <- match.arg(method)

  if (is.data.frame(x)) {
    if (!column %in% names(x)) {
      stop("`x` has no column named \"", column, "\".", call. = FALSE)
    }
    if (.matches_discrete_name(column, discrete_names)) {
      stop(
        "`column` (\"", column, "\") looks like discrete/categorical data ",
        "(matches a name in `discrete_names`); contour levels don't make sense ",
        "for it. Use `slice_label_contours()` instead, or pass a different `column`.",
        call. = FALSE
      )
    }
    values <- x[[column]]
  } else {
    if (is.character(x)) x <- ReadImage_fix(x)
    check_sitk_image(x)
    values <- as.vector(SimpleITK::as.array(x))
  }

  if (.looks_discrete(values)) {
    stop(
      "The sampled values look discrete/categorical (mostly integer-valued, ",
      "few unique values); contour levels don't make sense for them. ",
      "Use `slice_label_contours()` instead.",
      call. = FALSE
    )
  }

  if (method == "quantile") {
    if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n)) {
      stop("`n` must be a single positive integer.", call. = FALSE)
    }
    .quantile_levels(values, n)
  } else {
    if (!is.numeric(min_n) || length(min_n) != 1 || min_n < 1 || min_n != round(min_n)) {
      stop("`min_n` must be a single positive integer.", call. = FALSE)
    }
    .trough_levels(values, min_n)
  }
}
