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
