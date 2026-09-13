#' Flip a MINC image's x/y orientation for display, preserving pixel data
#'
#' @description
#' MINC files, when read via `SimpleITK`, often have direction cosines that
#' are the mirror image of what a NIfTI conversion of the same anatomy would
#' show. This flips the pixel data along x/y (via `SimpleITK::FlipImageFilter()`)
#' and recomputes the origin so the image occupies the same physical bounding
#' box, mirrored — a display-orientation fix, not a coordinate registration
#' between formats (it does not, and is not meant to, make a MINC-read image's
#' world coordinates numerically match an independently-converted NIfTI file
#' of the same anatomy; see the geometry-model notes in `CLAUDE.md` for the
#' real-data investigation confirming this).
#'
#' `z` is never flipped. Used by both [ReadImage_fix()] (to correct on read)
#' and [WriteImage_fix()] (to undo the correction on write, for MINC output).
#'
#' @param image A `SimpleITK` image.
#' @return The flipped `SimpleITK` image, with `OriginalFileType` metadata
#'   copied over if present.
#' @examples
#' \dontrun{
#' image <- SimpleITK::ReadImage("brain_image.mnc")
#' flipped <- orientation_correction(image)
#' }
#' @keywords internal
orientation_correction <- function(image) {
  flip_filter <- SimpleITK::FlipImageFilter()
  flip_filter$SetFlipAxes(c(TRUE, TRUE, FALSE)) # flip x and y, not z

  flipped_image <- flip_filter$Execute(image)

  if (image$HasMetaDataKey("OriginalFileType")) {
    flipped_image$SetMetaData("OriginalFileType", image$GetMetaData("OriginalFileType"))
  }

  flipped_image
}

#' Read an image file, correcting MINC orientation for display
#'
#' @description
#' Reads `file` via `SimpleITK::ReadImage()`. If `file` is a MINC file
#' (`.mnc`/`.minc`), also applies [orientation_correction()] and records the
#' original MINC direction matrix and file type as image metadata
#' (`OriginalFileType`, `OriginalDirection`), so [WriteImage_fix()] can later
#' undo the correction if writing back out to MINC or another format.
#'
#' @param file Path to an image file.
#' @return A `SimpleITK` image.
#' @examples
#' \dontrun{
#' # A MINC file gets its orientation corrected for display; a NIfTI file is
#' # returned unchanged (still tagged with metadata for WriteImage_fix()).
#' image <- ReadImage_fix("brain_image.mnc")
#' image2 <- ReadImage_fix("brain_image.nii")
#' }
#' @export
ReadImage_fix <- function(file) {
  image <- SimpleITK::ReadImage(file)
  is_minc <- tolower(tools::file_ext(file)) %in% c("mnc", "minc")

  if (!is_minc) {
    image$SetMetaData("OriginalFileType", "Other")
    return(image)
  }

  image$SetMetaData("OriginalFileType", "MINC")
  orig_direction <- image$GetDirection()
  image$SetMetaData("OriginalDirection", paste(orig_direction, collapse = ","))

  corrected_image <- orientation_correction(image)
  corrected_image$SetMetaData("OriginalFileType", "MINC")
  corrected_image$SetMetaData("OriginalDirection", paste(orig_direction, collapse = ","))

  corrected_image
}

#' Write an image file, restoring MINC orientation if needed
#'
#' @description
#' Writes `image` via `SimpleITK::WriteImage()`. If `image` was originally
#' read from a MINC file (tracked via the `OriginalFileType` metadata
#' [ReadImage_fix()] sets), the [orientation_correction()] applied on read is
#' undone before writing: fully (re-flipping) if the output is also MINC, or
#' by restoring the original MINC direction matrix (from the
#' `OriginalDirection` metadata) for any other output format. Images not
#' originally read as MINC are written unchanged.
#'
#' Note a real, pre-existing cross-language inconsistency, not addressed
#' here: the Python package's `WriteImage_fix()` only restores the original
#' direction for `.nii` output specifically (not `.nii.gz` or other non-MINC
#' formats), and hardcodes a generic RAS direction matrix rather than
#' restoring the image's actual original one. Left as-is pending a
#' deliberate decision to reconcile the two (a behavior change, not a pure
#' reorganization).
#'
#' @param image A `SimpleITK` image, typically from [ReadImage_fix()].
#' @param output_file Output file path.
#' @return `NULL`, invisibly; called for the side effect of writing `output_file`.
#' @examples
#' \dontrun{
#' image <- ReadImage_fix("brain_image.mnc")
#'
#' # MINC -> MINC: re-flipped back to proper MINC orientation before writing.
#' WriteImage_fix(image, "output.mnc")
#'
#' # MINC -> NIfTI: original MINC direction matrix restored before writing.
#' WriteImage_fix(image, "output.nii")
#' }
#' @export
WriteImage_fix <- function(image, output_file) {
  check_sitk_image(image)

  was_original_minc <- image$HasMetaDataKey("OriginalFileType") &&
    image$GetMetaData("OriginalFileType") == "MINC"

  output_extension <- tolower(tools::file_ext(output_file))
  is_output_minc <- output_extension %in% c("mnc", "minc")

  image_to_write <- image

  if (was_original_minc && is_output_minc) {
    # MINC -> MINC: undo the reading correction to restore proper MINC orientation.
    image_to_write <- orientation_correction(image)
  } else if (was_original_minc && !is_output_minc && image$HasMetaDataKey("OriginalDirection")) {
    # MINC -> any other format: restore the original MINC direction matrix directly.
    orig_dir <- as.numeric(strsplit(image$GetMetaData("OriginalDirection"), ",")[[1]])
    image_to_write$SetDirection(orig_dir)
  }

  SimpleITK::WriteImage(image_to_write, output_file)
  invisible(NULL)
}
