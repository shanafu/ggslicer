# Build a SimpleITK image of arbitrary dimension (2-5), with pixel values
# assigned by `value_fn(idx)`, where `idx` is the 0-indexed integer index
# vector for that pixel (length = length(size)).
filled_sitk_image <- function(size, value_fn, origin = NULL, spacing = NULL) {
  img <- SimpleITK::Image(as.integer(size), "sitkFloat32")
  if (!is.null(origin)) img$SetOrigin(origin)
  if (!is.null(spacing)) img$SetSpacing(spacing)

  axes <- lapply(size, function(n) 0:(n - 1))
  names(axes) <- paste0("d", seq_along(size))
  grid <- expand.grid(axes, KEEP.OUT.ATTRS = FALSE)

  for (r in seq_len(nrow(grid))) {
    idx <- as.integer(grid[r, ])
    img$SetPixel(idx, value_fn(idx))
  }
  img
}
