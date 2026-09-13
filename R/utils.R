#' Validate that `image` is a SimpleITK image object
#'
#' Gives a specifically helpful error if the caller passed a file path
#' instead (an easy mistake: several functions elsewhere in this package,
#' e.g. `slice_axis()`/`ReadImage_fix()`, take a path and read it
#' internally, but the functions that call this validator expect an
#' already-loaded image).
#'
#' @param image The value to check.
#' @param arg_name Name to use for `image` in the error message.
#'
#' @returns `TRUE`, invisibly, if `image` is a valid SimpleITK image; otherwise errors.
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

#' Internal function to read and fix orientation of MINC images
#'
#' @param image An image object.
#'
#' @returns Fixed image object.
#' @keywords internal
#'
#' @importFrom SimpleITK FlipImageFilter
orientation_correction <- function(image) {
  flip_filter <- FlipImageFilter()

  # Flip along axes 0 and 1 (first two dimensions)
  flip_axes <- c(T, T, F)  # Flip X and Y, not Z
  flip_filter$SetFlipAxes(flip_axes)

  # Pass image through flip_filter
  flipped_image <- flip_filter$Execute(image)

  # Copy metadata
  if (image$HasMetaDataKey("OriginalFileType")) {
    flipped_image$SetMetaData("OriginalFileType", image$GetMetaData("OriginalFileType"))
  }

  return(flipped_image)

}
