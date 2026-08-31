#' Read an image and when MINC image, flip image orientation while preserving all image data
#'
#' @param file Path to an image file.
#'
#' @returns An image object or corrected image object.
#' @export
#'
#' @examples
#'  \dontrun{
#' file <- "path/to/image"
#' ReadImage_fix(file)
#' }
#'
#' @importFrom SimpleITK ReadImage
#' @importFrom tools file_ext
ReadImage_fix <- function(file) {
  image <- ReadImage(file)
  extension <- tolower(file_ext(file))
  is_minc <- extension %in% c("mnc", "minc")

  if (is_minc) {
    image$SetMetaData("OriginalFileType", "MINC")

    # Store original MINC direction matrix
    orig_direction <- image$GetDirection()
    image$SetMetaData("OriginalDirection", paste(orig_direction, collapse=","))

    # Apply orientation correction for proper display
    corrected_image <- orientation_correction(image)

    # Copy metadata to corrected image
    corrected_image$SetMetaData("OriginalFileType", "MINC")
    corrected_image$SetMetaData("OriginalDirection", paste(orig_direction, collapse=","))

    image <- corrected_image

    return(image)

  } else {
    image$SetMetaData("OriginalFileType", "Other")
    return(image)
  }
}
