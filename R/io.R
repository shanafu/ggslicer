#' Convert a MINC image's header between MINC-native and RAS/LPS-consistent orientation
#'
#' @description
#' `SimpleITK::ReadImage()` handles NIfTI's stored orientation correctly: NIfTI
#' files conventionally store their qform/sform in a right-handed
#' Right-Anterior-Superior (RAS) convention, and SimpleITK converts this to
#' ITK's own internal Left-Posterior-Superior (LPS) convention on read (by
#' negating x/y), exactly as it does for DICOM. It does **not** perform the
#' equivalent conversion for MINC: a MINC file's own dimension metadata
#' (confirmed directly via `mincheader`/`mincinfo`, and by reading a MINC2
#' file's HDF5 structure directly, bypassing SimpleITK entirely) already
#' describes a RAS-like convention (e.g. `xspace`'s comment
#' "X increases from patient left to right", `yspace`'s
#' "Y increases from patient posterior to anterior", with positive
#' `direction_cosines` and a positive `step`) — but SimpleITK reports this
#' completely unconverted (direction reported as identity, origin exactly as
#' stored), leaving MINC-read images in a different, inconsistent coordinate
#' convention from NIfTI-read images of the very same anatomy.
#'
#' This function performs the missing conversion directly on the image's
#' header (origin and direction), by negating their `x`/`y` components (`z`
#' is never touched) — equivalent to left-multiplying both by
#' `diag(-1, -1, 1)`. **The pixel array is never reordered or copied in any
#' way that changes voxel content** — confirmed directly (comparing raw MINC
#' and NIfTI conversions of the same anatomy, for multiple real fixtures)
#' that the array is already in the same index order in both formats; only
#' the header's own description of that array's orientation was wrong for
#' MINC. This is a correction, not a data transformation, and it is its own
#' inverse: applying it twice returns the exact original header.
#'
#' (This function used to flip the pixel data via
#' `SimpleITK::FlipImageFilter()` and recompute the origin from the volume's
#' own bounding box. That was a real, confirmed bug: since the array never
#' needed reordering, physically flipping it produced an image whose header
#' was self-consistent but did not describe the same real-world locations as
#' an independently-converted NIfTI file of the same anatomy — verified
#' directly to disagree by tens of millimeters along `y` for a real fixture.
#' See `CLAUDE.md` for the full investigation, evidence, and fix.)
#'
#' @param image A `SimpleITK` image.
#' @return A new, independent `SimpleITK` image (the input is not modified)
#'   with `origin`/`direction` corrected, and `OriginalFileType` metadata
#'   copied over if present.
#' @examples
#' \dontrun{
#' image <- SimpleITK::ReadImage("brain_image.mnc")
#' corrected <- orientation_correction(image)
#' }
#' @keywords internal
orientation_correction <- function(image) {
  flip_diag <- diag(c(-1, -1, 1))

  origin <- image$GetOrigin()
  new_origin <- as.numeric(flip_diag %*% origin)

  # GetDirection()/SetDirection() flatten row-major; that matrix's columns
  # are each index axis's world direction (see CLAUDE.md's geometry-model
  # notes). Left-multiplying by flip_diag negates the x/y world-space
  # component of every column, i.e. of every index axis's direction vector.
  direction_mat <- matrix(image$GetDirection(), nrow = 3, byrow = TRUE)
  new_direction_mat <- flip_diag %*% direction_mat

  corrected <- SimpleITK::Image(image) # independent copy; pixel data is shared/COW, not reordered
  corrected$SetOrigin(new_origin)
  corrected$SetDirection(as.vector(t(new_direction_mat)))

  if (image$HasMetaDataKey("OriginalFileType")) {
    corrected$SetMetaData("OriginalFileType", image$GetMetaData("OriginalFileType"))
  }

  corrected
}

#' Read an image file, correcting MINC orientation for display
#'
#' @description
#' Reads `file` via `SimpleITK::ReadImage()`. If `file` is a MINC file
#' (`.mnc`/`.minc`), also applies [orientation_correction()] so that its
#' `origin`/`direction` describe the same real-world locations an
#' independently-converted NIfTI file of the same anatomy would (see that
#' function's docs for why this is needed at all). The corrected image is
#' tagged with `OriginalFileType` metadata (`"MINC"` or `"Other"`) so
#' [WriteImage_fix()] knows whether to undo the correction on write.
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
  corrected_image <- orientation_correction(image)
  corrected_image$SetMetaData("OriginalFileType", "MINC")
  corrected_image
}

#' Write an image file, restoring MINC orientation if needed
#'
#' @description
#' Writes `image` via `SimpleITK::WriteImage()`. If `image` was originally
#' read from a MINC file (tracked via the `OriginalFileType` metadata
#' [ReadImage_fix()] sets) and `output_file` is also MINC, the
#' [orientation_correction()] applied on read is undone (it is its own
#' inverse) before writing, so the output MINC file's header matches what a
#' native MINC reader expects. For any other output format, `image`'s
#' current (corrected) header is already the anatomically correct one to
#' write — no restoration is needed, unlike previous versions of this
#' function (see `CLAUDE.md`). Images not originally read as MINC are always
#' written unchanged.
#'
#' @param image A `SimpleITK` image, typically from [ReadImage_fix()].
#' @param output_file Output file path.
#' @return `NULL`, invisibly; called for the side effect of writing `output_file`.
#' @examples
#' \dontrun{
#' image <- ReadImage_fix("brain_image.mnc")
#'
#' # MINC -> MINC: re-corrected back to proper MINC orientation before writing.
#' WriteImage_fix(image, "output.mnc")
#'
#' # MINC -> NIfTI: already has the correct header; written as-is.
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
    # MINC -> MINC: undo the reading correction (self-inverse) to restore
    # proper native-MINC orientation.
    image_to_write <- orientation_correction(image)
  }

  SimpleITK::WriteImage(image_to_write, output_file)
  invisible(NULL)
}
