#' Extract voxel and world coordinates along specified slice
#'
#' @param file Path to an image file.
#' @param slice_orientation The desired axis to slice along. Options: "x"/"saggital"/"1", "y"/"coronal"/"2", "z"/"axial"/"3".
#' @param slice_orientation_coords Numeric vector of voxel coordinates specifying which slice planes to extract along the chosen axis.
#'
#' @returns A data frame with columns: x, y, z (world coordinates) and i, j, k (voxel coordinates) for all voxels on the specified slice planes.
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract coronal slices at y-coordinates 7.5 and 88
#' coords <- slice_axis("brain_image.nii", "y", c(7.5, 88))
#'
#' # Extract sagittal slice at x-coordinate 0
#' coords <- slice_axis("brain_image.nii", "sagittal", 0)
#'
#' # Extract axial slices
#' coords <- slice_axis("brain_image.nii", "z", c(-10, 0, 10))
#' }
#'
#' @importFrom pbapply pbapply
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange
#' @importFrom rlang sym
slice_axis <- function(file, slice_orientation, slice_orientation_coords) {

  image <- ReadImage_fix(file)
  file_dims <- image$GetSize()

  # Creates all possible voxel combinations with expand.grid
  voxel_coords <- expand.grid(
    i = 0:(file_dims[1] - 1),
    j = 0:(file_dims[2] - 1),
    k = 0:(file_dims[3] - 1)
  )

  # Convert voxel coordinates to world coordinates using the correct image
  world_coords_matrix <- pbapply(voxel_coords[, c("i", "j", "k")], 1, function(index) {
    image$TransformIndexToPhysicalPoint(as.numeric(index))
  })

  # Add world coordinates as columns
  voxel_coords$x <- world_coords_matrix[1,]
  voxel_coords$y <- world_coords_matrix[2,]
  voxel_coords$z <- world_coords_matrix[3,]

  axis_choice <- tolower(slice_orientation)
  filter_col_world <- NULL

  if (axis_choice %in% c("x", "sagittal", "1")) {
    filter_col_world <- "x"
  } else if (axis_choice %in% c("y", "coronal", "2")) {
    filter_col_world <- "y"
  } else if (axis_choice %in% c("z", "axial", "3")) {
    filter_col_world <- "z"
  } else {
    message("Invalid 'slice_orientation': '", slice_orientation, "'. Please choose 'x', 'sagittal', '1'; 'y', 'coronal', '2'; or 'z', 'axial', '3'.")
    return(NULL)
  }

  # Get unique world coordinates for the chosen axis
  available_world_coords <- unique(voxel_coords[[filter_col_world]])

  # Find closest available coordinates to user input
  target_slices <- sapply(slice_orientation_coords, function(target_coord) {
    available_world_coords[which.min(abs(available_world_coords - target_coord))]
  })

  target_slices <- unique(target_slices)

  # Loop through each target coordinate in the order the user provided them, filtering and appending to the running data frame
  tolerance <- 1e-6
  df_slice_axis <- NULL

  for (loop_index in seq_along(target_slices)) {

    target_coord <- target_slices[loop_index]

    df_temp <- voxel_coords %>%
      filter(abs(!!sym(filter_col_world) - target_coord) < tolerance) %>%
      mutate(slice_index = loop_index,
             slice_world = !!sym(filter_col_world),
             slice_axis = filter_col_world)

    if (is.null(df_slice_axis)) {
      df_slice_axis <- df_temp
    } else {
      df_slice_axis <- bind_rows(df_slice_axis, df_temp)
    }
  }

  # Reorder columns
  df_slice_axis <- df_slice_axis[, c("x", "y", "z", "i", "j", "k", "slice_axis", "slice_index", "slice_world")]

  return(df_slice_axis)
}
