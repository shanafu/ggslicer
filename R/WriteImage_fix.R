#' Write out an image file and when MINC image, flip image orientation while preserving all image data
#'
#' @param image An image object from ReadImage() or ReadImage_fix()
#' @param output_file Name of export file.
#'
#' @returns An exported image file.
#' @export
#'
#' @examples
#'  \dontrun{
#' # Write out MINC file
#' file <- "path/to/image"
#' image <- ReadImage_fix(file)
#' WriteImage_fix(image, "output.mnc")
#'
#' #Write out Nifti file
#' file <- "path/to/image"
#' image <- ReadImage_fix(file)
#' WriteImage_fix(image, "output.nii")
#' }
#'
#' @importFrom SimpleITK WriteImage
WriteImage_fix <- function(image, output_file) {
  # Check metadata
  original_file_type_meta <- image$HasMetaDataKey("OriginalFileType")
  was_original_minc <- F

  if (original_file_type_meta) {
    original_type <- image$GetMetaData("OriginalFileType")
    was_original_minc <- (original_type == "MINC")
  }

  output_extension <- tolower(file_ext(output_file))
  is_output_minc <- output_extension %in% c("mnc", "minc")

  image_to_write <- image

  if (was_original_minc && is_output_minc) {
    # MINC -> MINC: Undo the reading correction to restore proper MINC orientation
    image_to_write <- orientation_correction(image)

  } else if (was_original_minc && !is_output_minc) {
    # MINC -> non-MINC: Undo the reading correction to restore proper MINC orientation
    image_to_write <- image
    if (image$HasMetaDataKey("OriginalDirection")) {
      orig_dir_str <- image$GetMetaData("OriginalDirection")
      orig_dir <- as.numeric(strsplit(orig_dir_str, ",")[[1]])
      image_to_write$SetDirection(orig_dir)
    }
  }

  WriteImage(image_to_write, output_file)

}
