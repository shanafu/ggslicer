#' Append intensity values to slice_axis data frame
#'
#' @param file Path to an image file.
#' @param df_slice_axis Data frame containing voxel coordinates with columns: i, j, k (voxel indices) and x, y, z (world coordinates). Typically output from slice_axis() function. Default is NULL.
#'
#' @returns A data frame with columns: intensity, x, y, z, i, j, k.
#' @export
#'
#' @examples
#' \dontrun{
#' # First get coordinates from slice_axis
#' coords <- slice_axis("brain_image.nii", "y", c(7.5, 88))
#'
#' # Then extract intensity values at those coordinates
#' result <- slice_intensity("brain_image.nii", coords)}
slice_intensity <- function(file, df_slice_axis = NULL) {

  image <- ReadImage_fix(file) # Returns SimpleITK image

  # Get image dimensions and pizel data from SimpleITK image
  size <- image$GetSize() # Returns [width, height, depth]

  # Extract intensity values
  intensity_values <- numeric(nrow(df_slice_axis))

  for (row_index in 1:nrow(df_slice_axis)) {
    i_coord <- df_slice_axis$i[row_index]
    j_coord <- df_slice_axis$j[row_index]
    k_coord <- df_slice_axis$k[row_index]

    # Check bounds against image size
    if (i_coord >= 0 && i_coord < size[1] &&
        j_coord >= 0 && j_coord < size[2] &&
        k_coord >= 0 && k_coord < size[3]) {

      intensity_values[row_index] <- image$GetPixel(c(i_coord, j_coord, k_coord))
    } else {
      intensity_values[row_index] <- NA
    }
  }

  df_slice_axis$intensity <- intensity_values

  # Return properly ordered dataframe
  df <- df_slice_axis[, c("intensity", "x", "y", "z", "i", "j", "k", "slice_axis", "slice_index", "slice_world")]

  return(df)
}
