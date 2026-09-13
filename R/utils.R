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
