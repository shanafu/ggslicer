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
#'
#' @importFrom pbapply pbapply
slice_intensity <- function(file, df_slice_axis = NULL) {

  image <- ReadImage_fix(file) # Returns SimpleITK image

  # Get image dimensions and pizel data from SimpleITK image
  size <- image$GetSize() # Returns [width, height, depth]

  # Extract intensity values with progress bar
  intensity_values <- pbapply::pbapply(df_slice_axis[, c("i", "j", "k")], 1, function(row) {
    i_coord <- row[1]
    j_coord <- row[2]
    k_coord <- row[3]

    if (i_coord >= 0 && i_coord < size[1] &&
        j_coord >= 0 && j_coord < size[2] &&
        k_coord >= 0 && k_coord < size[3]) {
      image$GetPixel(as.numeric(c(i_coord, j_coord, k_coord)))
    } else {
      NA
    }
  })

  df_slice_axis$intensity <- intensity_values

  # Check if optional slice columns exist before selecting
  cols_to_select <- c("intensity", "x", "y", "z", "i", "j", "k", "slice_axis", "slice_index", "slice_world")
  available_cols <- intersect(cols_to_select, colnames(df_slice_axis))

  # Return properly ordered dataframe
  df <- df_slice_axis[, available_cols]
  
  return(df_slice_axis[, available_cols])
}
