# Internal numeric tolerance for unit-norm / orthogonality checks.
.slice_orthonormal_tol <- 1e-6

# Cross product of two length-3 vectors.
.cross3 <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

# Validate that direction_i/direction_j are unit vectors and orthogonal.
.validate_direction <- function(direction_i, direction_j, tol = .slice_orthonormal_tol) {
  if (!is.numeric(direction_i) || length(direction_i) != 3) {
    stop("`direction_i` must be a numeric length-3 vector.", call. = FALSE)
  }
  if (!is.numeric(direction_j) || length(direction_j) != 3) {
    stop("`direction_j` must be a numeric length-3 vector.", call. = FALSE)
  }
  if (abs(sqrt(sum(direction_i^2)) - 1) > tol) {
    stop("`direction_i` must be a unit vector (norm 1).", call. = FALSE)
  }
  if (abs(sqrt(sum(direction_j^2)) - 1) > tol) {
    stop("`direction_j` must be a unit vector (norm 1).", call. = FALSE)
  }
  if (abs(sum(direction_i * direction_j)) > tol) {
    stop("`direction_i` and `direction_j` must be orthogonal (dot product ~ 0).", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve an axis specification to an integer index 1/2/3. This is the
# standard axis-input convention for this package: numeric (1, 2, 3),
# Cartesian ("x", "y", "z"), or anatomical terms assuming right-anterior-
# superior (RAS+) positive orientation ("sagittal", "coronal", "axial", or
# "horizontal" as a synonym for axial). All axis-taking functions in this
# package should accept all three styles via this helper.
.resolve_axis_index <- function(axis) {
  if (is.numeric(axis)) {
    if (length(axis) != 1 || !(axis %in% c(1, 2, 3))) {
      stop("`axis` must be 1, 2, or 3 when given numerically.", call. = FALSE)
    }
    return(as.integer(axis))
  }
  if (!is.character(axis) || length(axis) != 1) {
    stop(
      "`axis` must be a single value: 1/2/3, \"x\"/\"y\"/\"z\", or an anatomical ",
      "term (\"sagittal\"/\"coronal\"/\"axial\"/\"horizontal\").",
      call. = FALSE
    )
  }
  choice <- tolower(axis)
  if (choice %in% c("x", "sagittal", "1")) {
    return(1L)
  }
  if (choice %in% c("y", "coronal", "2")) {
    return(2L)
  }
  if (choice %in% c("z", "axial", "horizontal", "3")) {
    return(3L)
  }
  stop(
    "Invalid `axis`; use 1/2/3, \"x\"/\"y\"/\"z\", or \"sagittal\"/\"coronal\"/\"axial\"/\"horizontal\".",
    call. = FALSE
  )
}

#' A rectangular, evenly-sampled 2D slice through 3D physical space
#'
#' @description
#' `SliceGeometry` represents the geometry of a rectangular 2D slice embedded in 3D
#' physical ("world") space: a starting point (`origin`), an orthonormal
#' in-plane basis (`direction_i`, `direction_j`), per-axis step sizes
#' (`spacing`), and per-axis sample counts (`size`). Sampling need not be
#' isotropic, and the slice need not be parallel to any Cartesian axis.
#'
#' `SliceGeometry` holds no pixel/intensity data itself and has no dependency on any
#' particular image — it is pure geometry, analogous to how DICOM describes a
#' single slice (`ImagePositionPatient`/`ImageOrientationPatient`/
#' `PixelSpacing`) or how a SimpleITK image stores its own
#' origin/direction/spacing/size.
#'
#' @export
SliceGeometry <- R6::R6Class(
  "SliceGeometry",
  public = list(

    #' @description Create a new `SliceGeometry`.
    #' @param origin Numeric length-3 vector: world coordinates of sample `(i = 0, j = 0)`.
    #' @param direction_i Numeric length-3 unit vector: world direction of the `i` axis.
    #' @param direction_j Numeric length-3 unit vector: world direction of the `j` axis. Must be orthogonal to `direction_i`.
    #' @param spacing Numeric length-2 vector: step size along `(i, j)`. Both entries must be `> 0`.
    #' @param size Integer-valued length-2 vector: number of samples along `(i, j)`. Both entries must be `>= 1`.
    initialize = function(origin, direction_i, direction_j, spacing, size) {
      private$set_origin_impl(origin)
      private$set_direction_impl(direction_i, direction_j)
      private$set_spacing_impl(spacing)
      private$set_size_impl(size)
    },

    #' @description World coordinates of sample `(i = 0, j = 0)`.
    get_origin = function() private$origin_,

    #' @description The 3x2 direction matrix (columns `i`, `j`; rows `x`, `y`, `z`).
    get_direction = function() private$direction_,

    #' @description World-space unit vector for the `i` axis.
    get_direction_i = function() private$direction_[, "i"],

    #' @description World-space unit vector for the `j` axis.
    get_direction_j = function() private$direction_[, "j"],

    #' @description Step size along `(i, j)`.
    get_spacing = function() private$spacing_,

    #' @description Number of samples along `(i, j)`.
    get_size = function() private$size_,

    #' @description Unit normal vector of the slice's plane (`direction_i x direction_j`).
    get_normal = function() .cross3(private$direction_[, "i"], private$direction_[, "j"]),

    #' @description The infinite plane the slice lies on, independent of its
    #'   finite extent: `list(point, normal)`.
    get_plane = function() list(point = private$origin_, normal = self$get_normal()),

    #' @description Physical extent (width, height) of the slice: `spacing * (size - 1)`.
    get_extent = function() private$spacing_ * (private$size_ - 1),

    #' @description World coordinates of the grid's midpoint.
    get_center = function() {
      extent <- self$get_extent()
      private$origin_ +
        0.5 * extent[1] * private$direction_[, "i"] +
        0.5 * extent[2] * private$direction_[, "j"]
    },

    #' @description World coordinates of the 4 corners of the sampling
    #'   rectangle, as a tibble.
    get_bounds = function() {
      extent <- self$get_extent()
      di <- private$direction_[, "i"]
      dj <- private$direction_[, "j"]
      corners <- rbind(
        private$origin_,
        private$origin_ + extent[1] * di,
        private$origin_ + extent[2] * dj,
        private$origin_ + extent[1] * di + extent[2] * dj
      )
      tibble::tibble(
        corner = c("i0_j0", "i1_j0", "i0_j1", "i1_j1"),
        x = corners[, 1], y = corners[, 2], z = corners[, 3]
      )
    },

    #' @description Physical coordinates of every sample point in the slice,
    #'   as a tibble with columns `i, j, x, y, z`. `i`/`j` are 0-indexed, in
    #'   the same convention used elsewhere in this package. The result is
    #'   memoized and recomputed only after a setter changes the geometry.
    get_sample_points = function() {
      if (!is.null(private$sample_points_cache_)) {
        return(private$sample_points_cache_)
      }

      n_i <- private$size_[1]
      n_j <- private$size_[2]
      ij <- expand.grid(i = 0:(n_i - 1), j = 0:(n_j - 1))

      scaled <- as.matrix(ij) %*% diag(private$spacing_, nrow = 2)
      world <- scaled %*% t(private$direction_)
      world <- sweep(world, 2, private$origin_, "+")

      out <- tibble::tibble(
        i = ij$i, j = ij$j,
        x = world[, 1], y = world[, 2], z = world[, 3]
      )
      private$sample_points_cache_ <- out
      out
    },

    #' @description Update the origin.
    #' @param origin Numeric length-3 vector.
    set_origin = function(origin) {
      private$set_origin_impl(origin)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the in-plane direction basis.
    #' @param direction_i Numeric length-3 unit vector, orthogonal to `direction_j`.
    #' @param direction_j Numeric length-3 unit vector, orthogonal to `direction_i`.
    set_direction = function(direction_i, direction_j) {
      private$set_direction_impl(direction_i, direction_j)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the per-axis step size.
    #' @param spacing Numeric length-2 vector, both entries `> 0`.
    set_spacing = function(spacing) {
      private$set_spacing_impl(spacing)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the per-axis sample count.
    #' @param size Integer-valued length-2 vector, both entries `>= 1`.
    set_size = function(size) {
      private$set_size_impl(size)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Reposition the rectangle so its midpoint is at `center`,
    #'   keeping direction, spacing, and size unchanged (updates `origin`).
    #' @param center Numeric length-3 vector.
    set_center = function(center) {
      if (!is.numeric(center) || length(center) != 3) {
        stop("`center` must be a numeric length-3 vector.", call. = FALSE)
      }
      extent <- self$get_extent()
      new_origin <- center -
        0.5 * extent[1] * private$direction_[, "i"] -
        0.5 * extent[2] * private$direction_[, "j"]
      private$set_origin_impl(new_origin)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the step size along `i` only.
    #' @param spacing_i Single number, `> 0`.
    set_spacing_i = function(spacing_i) {
      if (!is.numeric(spacing_i) || length(spacing_i) != 1 || spacing_i <= 0) {
        stop("`spacing_i` must be a single number > 0.", call. = FALSE)
      }
      private$set_spacing_impl(c(spacing_i, private$spacing_[2]))
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the step size along `j` only.
    #' @param spacing_j Single number, `> 0`.
    set_spacing_j = function(spacing_j) {
      if (!is.numeric(spacing_j) || length(spacing_j) != 1 || spacing_j <= 0) {
        stop("`spacing_j` must be a single number > 0.", call. = FALSE)
      }
      private$set_spacing_impl(c(private$spacing_[1], spacing_j))
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the sample count along `i` only.
    #' @param size_i Single integer, `>= 1`.
    set_size_i = function(size_i) {
      if (length(size_i) != 1 || size_i < 1 || size_i != round(size_i)) {
        stop("`size_i` must be a single integer >= 1.", call. = FALSE)
      }
      private$set_size_impl(c(size_i, private$size_[2]))
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the sample count along `j` only.
    #' @param size_j Single integer, `>= 1`.
    set_size_j = function(size_j) {
      if (length(size_j) != 1 || size_j < 1 || size_j != round(size_j)) {
        stop("`size_j` must be a single integer >= 1.", call. = FALSE)
      }
      private$set_size_impl(c(private$size_[1], size_j))
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Set the physical extent (width, height) directly, deriving
    #'   either `spacing` (keeping `size` fixed) or `size` (keeping `spacing`
    #'   fixed) to match.
    #' @param extent Numeric length-2 vector, both entries `>= 0`.
    #' @param adjust Which field to derive: `"spacing"` (default, changes
    #'   resolution to fit the new extent at the current sample count) or
    #'   `"size"` (changes sample count to fit the new extent at the current
    #'   resolution).
    set_extent = function(extent, adjust = c("spacing", "size")) {
      adjust <- match.arg(adjust)
      if (!is.numeric(extent) || length(extent) != 2 || any(extent < 0)) {
        stop("`extent` must be a non-negative numeric length-2 vector.", call. = FALSE)
      }
      if (adjust == "spacing") {
        if (any(private$size_ < 2)) {
          stop(
            "Cannot derive `spacing` from `extent` when `size` < 2 along an axis; ",
            "use `adjust = \"size\"` or call `set_spacing()` directly.",
            call. = FALSE
          )
        }
        private$set_spacing_impl(extent / (private$size_ - 1))
      } else {
        private$set_size_impl(round(extent / private$spacing_) + 1)
      }
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Return a new `SliceGeometry`, identical to this one but with its
    #'   origin shifted by `offset` (e.g. `offset = spacing_k * get_normal()`
    #'   to move to the next parallel slice in a stack).
    #' @param offset Numeric length-3 vector.
    translate = function(offset) {
      if (!is.numeric(offset) || length(offset) != 3) {
        stop("`offset` must be a numeric length-3 vector.", call. = FALSE)
      }
      SliceGeometry$new(
        origin = private$origin_ + offset,
        direction_i = private$direction_[, "i"],
        direction_j = private$direction_[, "j"],
        spacing = private$spacing_,
        size = private$size_
      )
    },

    #' @description Build a pixel-less, single-voxel-thick SimpleITK
    #'   reference image whose geometry (origin/direction/spacing/size)
    #'   exactly matches this slice, suitable as the `referenceImage`
    #'   argument to `SimpleITK::Resample()` for extracting intensities.
    #' @param spacing_k Spacing to assign to the synthetic third (normal) axis. Arbitrary, since that axis has only one sample; defaults to `1`.
    as_sitk_reference_image = function(spacing_k = 1) {
      direction_3x3 <- cbind(private$direction_, normal = self$get_normal())
      img <- SimpleITK::Image(
        as.integer(private$size_[1]), as.integer(private$size_[2]), 1L,
        "sitkFloat32"
      )
      img$SetOrigin(private$origin_)
      img$SetSpacing(c(private$spacing_, spacing_k))
      img$SetDirection(as.vector(t(direction_3x3)))
      img
    },

    #' @description Print a short summary of the slice's geometry.
    #' @param ... Unused; present for compatibility with the generic `print()`.
    print = function(...) {
      cat("<SliceGeometry>\n")
      cat("  origin:     ", paste(signif(private$origin_, 4), collapse = ", "), "\n")
      cat("  direction_i:", paste(signif(private$direction_[, "i"], 4), collapse = ", "), "\n")
      cat("  direction_j:", paste(signif(private$direction_[, "j"], 4), collapse = ", "), "\n")
      cat("  normal:     ", paste(signif(self$get_normal(), 4), collapse = ", "), "\n")
      cat("  spacing:    ", paste(private$spacing_, collapse = ", "), "\n")
      cat("  size:       ", paste(private$size_, collapse = ", "), "\n")
      cat("  extent:     ", paste(signif(self$get_extent(), 4), collapse = ", "), "\n")
      invisible(self)
    }
  ),
  private = list(
    origin_ = NULL,
    direction_ = NULL,
    spacing_ = NULL,
    size_ = NULL,
    sample_points_cache_ = NULL,

    set_origin_impl = function(origin) {
      if (!is.numeric(origin) || length(origin) != 3) {
        stop("`origin` must be a numeric length-3 vector.", call. = FALSE)
      }
      private$origin_ <- as.numeric(origin)
    },

    set_direction_impl = function(direction_i, direction_j) {
      .validate_direction(direction_i, direction_j)
      direction <- cbind(as.numeric(direction_i), as.numeric(direction_j))
      dimnames(direction) <- list(c("x", "y", "z"), c("i", "j"))
      private$direction_ <- direction
    },

    set_spacing_impl = function(spacing) {
      if (!is.numeric(spacing) || length(spacing) != 2 || any(spacing <= 0)) {
        stop("`spacing` must be a numeric length-2 vector with both entries > 0.", call. = FALSE)
      }
      private$spacing_ <- as.numeric(spacing)
    },

    set_size_impl = function(size) {
      if (length(size) != 2 || any(size < 1) || any(size != round(size))) {
        stop("`size` must be an integer-valued length-2 vector with both entries >= 1.", call. = FALSE)
      }
      private$size_ <- as.integer(size)
    },

    invalidate_cache = function() {
      private$sample_points_cache_ <- NULL
    }
  )
)

#' @description
#' Construct a `SliceGeometry` from two opposite corners of the sampling rectangle.
#'
#' Two opposite corners alone do not uniquely determine a rectangle in 3D:
#' for a fixed diagonal `p1 - p0`, any orthogonal decomposition of that
#' diagonal into two edge vectors gives a different, equally valid rectangle
#' (differing in orientation and, often, in which plane it lies on). Supplying
#' `direction_i` resolves this: the diagonal is projected onto `direction_i`
#' to get the extent along `i`, and `direction_j` is derived as the
#' (automatically orthogonal) remainder, which also fixes the plane.
#'
#' @param p0 Numeric length-3 vector: world coordinates of the corner at `(i = 0, j = 0)`.
#' @param p1 Numeric length-3 vector: world coordinates of the opposite corner, at `(i = max, j = max)`.
#' @param direction_i Numeric length-3 unit vector, pointing from `p0` toward `p1` along the `i` axis.
#' @param spacing Numeric length-2 vector: step size along `(i, j)`. Exactly one of `spacing`/`size` must be supplied.
#' @param size Integer-valued length-2 vector: number of samples along `(i, j)`. Exactly one of `spacing`/`size` must be supplied.
#' @return A new `SliceGeometry`.
#' @rdname SliceGeometry
#' @name SliceGeometry_from_corners
SliceGeometry$from_corners <- function(p0, p1, direction_i, spacing = NULL, size = NULL) {
  if (is.null(spacing) == is.null(size)) {
    stop("Specify exactly one of `spacing` or `size`.", call. = FALSE)
  }
  if (!is.numeric(p0) || length(p0) != 3 || !is.numeric(p1) || length(p1) != 3) {
    stop("`p0` and `p1` must be numeric length-3 vectors.", call. = FALSE)
  }
  if (!is.numeric(direction_i) || length(direction_i) != 3) {
    stop("`direction_i` must be a numeric length-3 vector.", call. = FALSE)
  }
  if (abs(sqrt(sum(direction_i^2)) - 1) > .slice_orthonormal_tol) {
    stop("`direction_i` must be a unit vector (norm 1).", call. = FALSE)
  }

  d <- as.numeric(p1) - as.numeric(p0)
  w <- sum(d * direction_i)
  if (w <= 0) {
    stop(
      "`direction_i` must point from `p0` toward `p1` (projected extent must be positive).",
      call. = FALSE
    )
  }
  r <- d - w * direction_i
  h <- sqrt(sum(r^2))
  if (h <= .slice_orthonormal_tol) {
    stop(
      "`p0`, `p1`, and `direction_i` are collinear; no unique rectangle exists.",
      call. = FALSE
    )
  }
  direction_j <- r / h

  if (!is.null(size)) {
    if (length(size) != 2) stop("`size` must be a length-2 vector.", call. = FALSE)
    if (any(size < 2)) {
      stop(
        "`size` must be >= 2 along both axes to derive `spacing` from the corners; ",
        "supply `spacing` explicitly for a single-sample axis.",
        call. = FALSE
      )
    }
    spacing <- c(w, h) / (size - 1)
  } else {
    if (length(spacing) != 2) stop("`spacing` must be a length-2 vector.", call. = FALSE)
    size <- round(c(w, h) / spacing) + 1
  }

  SliceGeometry$new(
    origin = p0,
    direction_i = direction_i,
    direction_j = direction_j,
    spacing = spacing,
    size = size
  )
}

#' @description
#' Construct a `SliceGeometry` from its center point rather than its
#' `(i = 0, j = 0)` corner.
#'
#' @param center Numeric length-3 vector: world coordinates of the rectangle's midpoint.
#' @param direction_i Numeric length-3 unit vector: world direction of the `i` axis.
#' @param direction_j Numeric length-3 unit vector: world direction of the `j` axis. Must be orthogonal to `direction_i`.
#' @param spacing Numeric length-2 vector: step size along `(i, j)`.
#' @param size Integer-valued length-2 vector: number of samples along `(i, j)`.
#' @return A new `SliceGeometry`.
#' @rdname SliceGeometry
#' @name SliceGeometry_from_center
SliceGeometry$from_center <- function(center, direction_i, direction_j, spacing, size) {
  s <- SliceGeometry$new(
    origin = c(0, 0, 0), direction_i = direction_i, direction_j = direction_j,
    spacing = spacing, size = size
  )
  s$set_center(center)
  s
}

#' @description
#' Construct a `SliceGeometry` from a plane normal rather than an explicit
#' in-plane basis.
#'
#' A normal alone does not fix the in-plane rotation (any rotation of
#' `direction_i`/`direction_j` about the normal is still a valid orthonormal
#' basis for the same plane). If `direction_i` is not supplied, a default seed
#' vector (the world x-axis, or the y-axis if the normal is nearly parallel to
#' x) is projected into the plane to pick one. If `direction_i` is supplied,
#' it is projected into the plane the same way, so it need not already be
#' exactly orthogonal to `normal`.
#'
#' @param origin Numeric length-3 vector: world coordinates of sample `(i = 0, j = 0)`.
#' @param normal Numeric length-3 vector, nonzero (need not be unit length).
#' @param spacing Numeric length-2 vector: step size along `(i, j)`.
#' @param size Integer-valued length-2 vector: number of samples along `(i, j)`.
#' @param direction_i Optional numeric length-3 vector used to seed the
#'   in-plane rotation (see Details); must not be parallel to `normal`.
#'   Defaults to a world-axis seed.
#' @return A new `SliceGeometry`.
#' @rdname SliceGeometry
#' @name SliceGeometry_from_normal
SliceGeometry$from_normal <- function(origin, normal, spacing, size, direction_i = NULL) {
  if (!is.numeric(normal) || length(normal) != 3) {
    stop("`normal` must be a numeric length-3 vector.", call. = FALSE)
  }
  normal_norm <- sqrt(sum(normal^2))
  if (normal_norm <= .slice_orthonormal_tol) {
    stop("`normal` must be nonzero.", call. = FALSE)
  }
  normal_unit <- normal / normal_norm

  if (is.null(direction_i)) {
    seed <- c(1, 0, 0)
    if (abs(sum(seed * normal_unit)) > 1 - .slice_orthonormal_tol) seed <- c(0, 1, 0)
  } else {
    if (!is.numeric(direction_i) || length(direction_i) != 3) {
      stop("`direction_i` must be a numeric length-3 vector.", call. = FALSE)
    }
    seed <- direction_i
  }

  proj <- seed - sum(seed * normal_unit) * normal_unit
  proj_norm <- sqrt(sum(proj^2))
  if (proj_norm <= .slice_orthonormal_tol) {
    stop("`direction_i` is parallel to `normal`; supply a different seed direction.", call. = FALSE)
  }
  di <- proj / proj_norm
  dj <- .cross3(normal_unit, di)

  SliceGeometry$new(origin = origin, direction_i = di, direction_j = dj, spacing = spacing, size = size)
}

#' @description
#' Construct a `SliceGeometry` matching an axis-aligned slice through an
#' existing 3D SimpleITK image — the geometric bridge to this package's
#' original `slice_axis()` semantics, generalized to an oblique image
#' direction matrix.
#'
#' The in-plane axes, their spacing, and their sample counts are taken
#' directly from `image`'s own geometry (the two axes other than `axis`); the
#' out-of-plane position is snapped to the nearest voxel plane to
#' `coordinate` along `axis`, exactly as the legacy `slice_axis()` did.
#'
#' @param image A 3D `SimpleITK` image.
#' @param axis Which image axis is out-of-plane: `1`/`2`/`3`, `"x"`/`"y"`/`"z"`,
#'   or (assuming right-anterior-superior orientation)
#'   `"sagittal"`/`"coronal"`/`"axial"`/`"horizontal"`.
#' @param coordinate Single number: the desired world coordinate along `axis`
#'   (snapped to the nearest voxel plane).
#' @return A new `SliceGeometry`.
#' @rdname SliceGeometry
#' @name SliceGeometry_from_image_axis
SliceGeometry$from_image_axis <- function(image, axis, coordinate) {
  if (image$GetDimension() != 3) {
    stop("`image` must be a 3D image.", call. = FALSE)
  }
  axis_index <- .resolve_axis_index(axis)

  size_full <- image$GetSize()
  spacing_full <- image$GetSpacing()
  direction_mat <- matrix(image$GetDirection(), nrow = 3, byrow = TRUE)
  in_plane <- setdiff(1:3, axis_index)

  n_along <- size_full[axis_index]
  world_along <- vapply(0:(n_along - 1), function(v) {
    idx <- c(0L, 0L, 0L)
    idx[axis_index] <- v
    image$TransformIndexToPhysicalPoint(as.integer(idx))[axis_index]
  }, numeric(1))
  nearest_v <- (0:(n_along - 1))[which.min(abs(world_along - coordinate))]
  fixed_idx <- c(0L, 0L, 0L)
  fixed_idx[axis_index] <- nearest_v
  origin <- image$TransformIndexToPhysicalPoint(as.integer(fixed_idx))

  SliceGeometry$new(
    origin = origin,
    direction_i = direction_mat[, in_plane[1]],
    direction_j = direction_mat[, in_plane[2]],
    spacing = c(spacing_full[in_plane[1]], spacing_full[in_plane[2]]),
    size = c(size_full[in_plane[1]], size_full[in_plane[2]])
  )
}

#' @description
#' Construct a `SliceGeometry` from the 4 labeled corners produced by
#' `get_bounds()` — the inverse of that getter. Unlike `from_corners()`, no
#' extra disambiguating direction is needed: 4 labeled corners fully
#' determine the rectangle (`direction_i`/`direction_j` come directly from
#' two of its edges).
#'
#' @param bounds A data frame/tibble with columns `corner, x, y, z`,
#'   containing rows labeled `"i0_j0"`, `"i1_j0"`, and `"i0_j1"` (as produced
#'   by `get_bounds()`). If a `"i1_j1"` row is also present, it is checked for
#'   consistency with the other three.
#' @param spacing Numeric length-2 vector: step size along `(i, j)`. Exactly one of `spacing`/`size` must be supplied.
#' @param size Integer-valued length-2 vector: number of samples along `(i, j)`. Exactly one of `spacing`/`size` must be supplied.
#' @return A new `SliceGeometry`.
#' @rdname SliceGeometry
#' @name SliceGeometry_from_bounds
SliceGeometry$from_bounds <- function(bounds, spacing = NULL, size = NULL) {
  if (is.null(spacing) == is.null(size)) {
    stop("Specify exactly one of `spacing` or `size`.", call. = FALSE)
  }
  required <- c("i0_j0", "i1_j0", "i0_j1")
  if (is.null(bounds$corner) || !all(required %in% bounds$corner)) {
    stop(
      "`bounds` must contain rows labeled \"i0_j0\", \"i1_j0\", and \"i0_j1\" ",
      "(as produced by `get_bounds()`).",
      call. = FALSE
    )
  }
  corner_xyz <- function(label) {
    row <- bounds[bounds$corner == label, c("x", "y", "z"), drop = FALSE]
    if (nrow(row) != 1) {
      stop("`bounds` must have exactly one row for corner \"", label, "\".", call. = FALSE)
    }
    as.numeric(unlist(row[1, ]))
  }
  p00 <- corner_xyz("i0_j0")
  p10 <- corner_xyz("i1_j0")
  p01 <- corner_xyz("i0_j1")

  edge_i <- p10 - p00
  edge_j <- p01 - p00
  extent_i <- sqrt(sum(edge_i^2))
  extent_j <- sqrt(sum(edge_j^2))
  if (extent_i <= .slice_orthonormal_tol || extent_j <= .slice_orthonormal_tol) {
    stop("Degenerate `bounds`: corners \"i0_j0\"/\"i1_j0\"/\"i0_j1\" must not coincide.", call. = FALSE)
  }
  direction_i <- edge_i / extent_i
  direction_j <- edge_j / extent_j

  if ("i1_j1" %in% bounds$corner) {
    p11 <- corner_xyz("i1_j1")
    expected_p11 <- p00 + edge_i + edge_j
    if (max(abs(p11 - expected_p11)) > 1e-4 * max(1, extent_i, extent_j)) {
      stop(
        "Corner \"i1_j1\" is inconsistent with \"i0_j0\"/\"i1_j0\"/\"i0_j1\"; ",
        "`bounds` does not describe a rectangle.",
        call. = FALSE
      )
    }
  }

  if (!is.null(size)) {
    if (length(size) != 2) stop("`size` must be a length-2 vector.", call. = FALSE)
    if (any(size < 2)) {
      stop(
        "`size` must be >= 2 along both axes to derive `spacing` from `bounds`; ",
        "supply `spacing` explicitly for a single-sample axis.",
        call. = FALSE
      )
    }
    spacing <- c(extent_i, extent_j) / (size - 1)
  } else {
    if (length(spacing) != 2) stop("`spacing` must be a length-2 vector.", call. = FALSE)
    size <- round(c(extent_i, extent_j) / spacing) + 1
  }

  SliceGeometry$new(
    origin = p00, direction_i = direction_i, direction_j = direction_j,
    spacing = spacing, size = size
  )
}

#' A stack of parallel, evenly-spaced 2D slices
#'
#' @description
#' `SlicePackage` represents a set of parallel [SliceGeometry] rectangles stacked
#' along their shared normal direction at a fixed spacing — a regular
#' rectangular-cuboid sampling volume, exactly analogous to how a stack of 2D
#' DICOM slices forms a 3D volume. It is built from a single base `SliceGeometry`
#' (the slice at `k = 0`) plus a step size and sample count along the normal.
#'
#' @export
SlicePackage <- R6::R6Class(
  "SlicePackage",
  public = list(

    #' @description Create a new `SlicePackage`.
    #' @param base_slice A [SliceGeometry] object: the slice at `k = 0`.
    #' @param spacing_k Step size along the normal direction. Must be `> 0`.
    #' @param size_k Number of parallel slices in the stack. Must be `>= 1`.
    initialize = function(base_slice, spacing_k, size_k) {
      private$set_base_slice_impl(base_slice)
      private$set_spacing_k_impl(spacing_k)
      private$set_size_k_impl(size_k)
    },

    #' @description A clone of the `SliceGeometry` at `k = 0`. A clone (not
    #'   the live internal object) is returned so that mutating it can't
    #'   silently desynchronize this package's cached sample points; use
    #'   `set_base_slice()` to actually change it.
    get_base_slice = function() private$base_slice_$clone(),

    #' @description Step size along the normal (`k`) direction.
    get_spacing_k = function() private$spacing_k_,

    #' @description Number of parallel slices in the stack.
    get_size_k = function() private$size_k_,

    #' @description Full 3D spacing, `c(spacing_i, spacing_j, spacing_k)`.
    get_spacing = function() c(private$base_slice_$get_spacing(), private$spacing_k_),

    #' @description Full 3D size, `c(n_i, n_j, n_k)`.
    get_size = function() c(private$base_slice_$get_size(), private$size_k_),

    #' @description Unit normal of the base slice's plane (shared by every slice in the stack).
    get_normal = function() private$base_slice_$get_normal(),

    #' @description The `k`-th `SliceGeometry` in the stack.
    #' @param k Single integer in `0:(size_k - 1)`.
    get_slice = function(k) {
      if (length(k) != 1 || k != round(k) || k < 0 || k >= private$size_k_) {
        stop("`k` must be a single integer in 0:(size_k - 1).", call. = FALSE)
      }
      private$base_slice_$translate(k * private$spacing_k_ * self$get_normal())
    },

    #' @description Physical coordinates of every sample point in the stack,
    #'   as a tibble with columns `i, j, k, x, y, z`. Memoized like
    #'   [SliceGeometry]'s `get_sample_points()`.
    get_sample_points = function() {
      if (!is.null(private$sample_points_cache_)) {
        return(private$sample_points_cache_)
      }

      base <- private$base_slice_
      base_size <- base$get_size()
      n_i <- base_size[1]
      n_j <- base_size[2]
      n_k <- private$size_k_
      ijk <- expand.grid(i = 0:(n_i - 1), j = 0:(n_j - 1), k = 0:(n_k - 1))

      spacing3 <- c(base$get_spacing(), private$spacing_k_)
      direction3 <- cbind(base$get_direction(), normal = self$get_normal())

      scaled <- as.matrix(ijk) %*% diag(spacing3, nrow = 3)
      world <- scaled %*% t(direction3)
      world <- sweep(world, 2, base$get_origin(), "+")

      out <- tibble::tibble(
        i = ijk$i, j = ijk$j, k = ijk$k,
        x = world[, 1], y = world[, 2], z = world[, 3]
      )
      private$sample_points_cache_ <- out
      out
    },

    #' @description Build a pixel-less 3D SimpleITK reference image spanning
    #'   the whole stack in a single geometry, suitable as the
    #'   `referenceImage` argument to `SimpleITK::Resample()`.
    as_sitk_reference_image = function() {
      base <- private$base_slice_
      size <- self$get_size()
      spacing <- self$get_spacing()
      direction3x3 <- cbind(base$get_direction(), normal = self$get_normal())

      img <- SimpleITK::Image(
        as.integer(size[1]), as.integer(size[2]), as.integer(size[3]),
        "sitkFloat32"
      )
      img$SetOrigin(base$get_origin())
      img$SetSpacing(spacing)
      img$SetDirection(as.vector(t(direction3x3)))
      img
    },

    #' @description Resample a 3D SimpleITK image onto this stack's geometry
    #'   in a single `SimpleITK::Resample()` call, returning
    #'   `get_sample_points()` with an added `intensity` column (`NA` outside
    #'   the source image's bounds).
    #' @param image A 3D `SimpleITK` image.
    #' @param interpolator SimpleITK interpolator name, e.g. `"sitkLinear"` (default) or `"sitkNearestNeighbor"`.
    sample_intensity = function(image, interpolator = "sitkLinear") {
      if (image$GetDimension() != 3) {
        stop(
          "`image` must be a 3D image; extract a spatial sub-volume first for ",
          "higher-dimensional images (this is what `SlicePackageSet$sample_intensity()` ",
          "does automatically).",
          call. = FALSE
        )
      }
      ref <- self$as_sitk_reference_image()
      resampled <- SimpleITK::Resample(image, ref, SimpleITK::Transform(), interpolator, NaN)
      intensity <- as.vector(SimpleITK::as.array(resampled))
      intensity[is.nan(intensity)] <- NA_real_

      out <- self$get_sample_points()
      out$intensity <- intensity
      out
    },

    #' @description Update the step size along the normal.
    #' @param spacing_k Single number, `> 0`.
    set_spacing_k = function(spacing_k) {
      private$set_spacing_k_impl(spacing_k)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the number of parallel slices.
    #' @param size_k Single integer, `>= 1`.
    set_size_k = function(size_k) {
      private$set_size_k_impl(size_k)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Replace the base (`k = 0`) slice.
    #' @param base_slice A [SliceGeometry] object.
    set_base_slice = function(base_slice) {
      private$set_base_slice_impl(base_slice)
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the full 3D spacing at once.
    #' @param spacing Numeric length-3 vector, `c(spacing_i, spacing_j, spacing_k)`.
    set_spacing = function(spacing) {
      if (!is.numeric(spacing) || length(spacing) != 3) {
        stop("`spacing` must be a numeric length-3 vector.", call. = FALSE)
      }
      private$base_slice_$set_spacing(spacing[1:2])
      private$set_spacing_k_impl(spacing[3])
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Update the full 3D size at once.
    #' @param size Integer-valued length-3 vector, `c(n_i, n_j, n_k)`.
    set_size = function(size) {
      if (length(size) != 3) {
        stop("`size` must be an integer-valued length-3 vector.", call. = FALSE)
      }
      private$base_slice_$set_size(size[1:2])
      private$set_size_k_impl(size[3])
      private$invalidate_cache()
      invisible(self)
    },

    #' @description Print a short summary of the package's geometry.
    #' @param ... Unused; present for compatibility with the generic `print()`.
    print = function(...) {
      cat("<SlicePackage>\n")
      cat("  size:        ", paste(self$get_size(), collapse = ", "), "\n")
      cat("  spacing:     ", paste(self$get_spacing(), collapse = ", "), "\n")
      cat("  normal:      ", paste(signif(self$get_normal(), 4), collapse = ", "), "\n")
      cat("  base origin: ", paste(signif(private$base_slice_$get_origin(), 4), collapse = ", "), "\n")
      cat("  direction_i: ", paste(signif(private$base_slice_$get_direction_i(), 4), collapse = ", "), "\n")
      cat("  direction_j: ", paste(signif(private$base_slice_$get_direction_j(), 4), collapse = ", "), "\n")
      invisible(self)
    }
  ),
  private = list(
    base_slice_ = NULL,
    spacing_k_ = NULL,
    size_k_ = NULL,
    sample_points_cache_ = NULL,

    set_base_slice_impl = function(base_slice) {
      if (!inherits(base_slice, "SliceGeometry")) {
        stop("`base_slice` must be a `SliceGeometry` object.", call. = FALSE)
      }
      private$base_slice_ <- base_slice
    },
    set_spacing_k_impl = function(spacing_k) {
      if (!is.numeric(spacing_k) || length(spacing_k) != 1 || spacing_k <= 0) {
        stop("`spacing_k` must be a single number > 0.", call. = FALSE)
      }
      private$spacing_k_ <- as.numeric(spacing_k)
    },
    set_size_k_impl = function(size_k) {
      if (length(size_k) != 1 || size_k < 1 || size_k != round(size_k)) {
        stop("`size_k` must be a single integer >= 1.", call. = FALSE)
      }
      private$size_k_ <- as.integer(size_k)
    },
    invalidate_cache = function() {
      private$sample_points_cache_ <- NULL
    }
  )
)

#' @description
#' Construct a `SlicePackage` from a total stack thickness rather than a
#' per-step spacing.
#'
#' @param base_slice A [SliceGeometry] object: the slice at `k = 0`.
#' @param extent_k Total physical thickness of the stack (from the first to
#'   the last slice). `spacing_k` is derived as `extent_k / (size_k - 1)`.
#' @param size_k Number of parallel slices in the stack. Must be `>= 2`.
#' @return A new `SlicePackage`.
#' @rdname SlicePackage
#' @name SlicePackage_from_extent_k
SlicePackage$from_extent_k <- function(base_slice, extent_k, size_k) {
  if (length(size_k) != 1 || size_k < 2 || size_k != round(size_k)) {
    stop(
      "`size_k` must be a single integer >= 2 to derive `spacing_k` from `extent_k`; ",
      "use the primary constructor directly for size_k = 1.",
      call. = FALSE
    )
  }
  if (!is.numeric(extent_k) || length(extent_k) != 1 || extent_k <= 0) {
    stop("`extent_k` must be a single number > 0.", call. = FALSE)
  }
  SlicePackage$new(base_slice = base_slice, spacing_k = extent_k / (size_k - 1), size_k = size_k)
}

#' @description
#' Construct a `SlicePackage` treating the given slice as the *middle* of the
#' stack, rather than as its first (`k = 0`) slice.
#'
#' @param center_slice A [SliceGeometry] object: the slice at the middle of the stack.
#' @param spacing_k Step size along the normal direction. Must be `> 0`.
#' @param size_k Number of parallel slices in the stack. Must be `>= 1`.
#' @return A new `SlicePackage`.
#' @rdname SlicePackage
#' @name SlicePackage_from_center_k
SlicePackage$from_center_k <- function(center_slice, spacing_k, size_k) {
  if (!is.numeric(spacing_k) || length(spacing_k) != 1 || spacing_k <= 0) {
    stop("`spacing_k` must be a single number > 0.", call. = FALSE)
  }
  if (length(size_k) != 1 || size_k < 1 || size_k != round(size_k)) {
    stop("`size_k` must be a single integer >= 1.", call. = FALSE)
  }
  half_stack_extent <- spacing_k * (size_k - 1) / 2
  base <- center_slice$translate(-half_stack_extent * center_slice$get_normal())
  SlicePackage$new(base_slice = base, spacing_k = spacing_k, size_k = size_k)
}

#' @description
#' Construct a `SlicePackage` by adopting a list of already-built,
#' individually-defined [SliceGeometry] objects (e.g. one per DICOM slice),
#' rather than generating a regular stack from a single base slice and a
#' step size. Validates that every slice shares the same
#' direction/spacing/size and that consecutive slices (in the order given)
#' are evenly spaced along their shared normal.
#'
#' @param slices A list of 2 or more [SliceGeometry] objects, ordered from `k = 0` onward.
#' @return A new `SlicePackage`.
#' @rdname SlicePackage
#' @name SlicePackage_from_slices
SlicePackage$from_slices <- function(slices) {
  if (!is.list(slices) || length(slices) < 2) {
    stop("`slices` must be a list of at least 2 `SliceGeometry` objects.", call. = FALSE)
  }
  if (!all(vapply(slices, inherits, logical(1), what = "SliceGeometry"))) {
    stop("Every element of `slices` must be a `SliceGeometry` object.", call. = FALSE)
  }

  base <- slices[[1]]
  ref_direction <- base$get_direction()
  ref_spacing <- base$get_spacing()
  ref_size <- base$get_size()

  for (i in seq_along(slices)[-1]) {
    s <- slices[[i]]
    if (!isTRUE(all.equal(s$get_direction(), ref_direction, check.attributes = FALSE))) {
      stop("All slices must share the same direction (be parallel); slice ", i, " does not.", call. = FALSE)
    }
    if (!isTRUE(all.equal(s$get_spacing(), ref_spacing))) {
      stop("All slices must share the same spacing; slice ", i, " does not.", call. = FALSE)
    }
    if (!isTRUE(all.equal(s$get_size(), ref_size))) {
      stop("All slices must share the same size; slice ", i, " does not.", call. = FALSE)
    }
  }

  normal <- unname(base$get_normal())
  base_origin <- base$get_origin()
  offsets <- vapply(slices, function(s) sum((s$get_origin() - base_origin) * normal), numeric(1))

  for (i in seq_along(slices)) {
    expected_origin <- base_origin + offsets[i] * normal
    if (max(abs(slices[[i]]$get_origin() - expected_origin)) > 1e-4) {
      stop("Slice ", i, "'s origin is not aligned along the shared normal direction.", call. = FALSE)
    }
  }
  step_sizes <- diff(offsets)
  spacing_k <- step_sizes[1]
  if (spacing_k <= 0 || !isTRUE(all.equal(step_sizes, rep(spacing_k, length(step_sizes)), tolerance = 1e-4))) {
    stop(
      "`slices` must be evenly spaced, in order, along their shared normal direction.",
      call. = FALSE
    )
  }

  SlicePackage$new(base_slice = base, spacing_k = spacing_k, size_k = length(slices))
}

#' @description
#' Construct a `SlicePackage` spanning an entire 3D SimpleITK image along one
#' of its axes, using the image's own resolution along that axis — the
#' `SlicePackage` equivalent of `SliceGeometry$from_image_axis()`.
#'
#' @param image A 3D `SimpleITK` image.
#' @param axis Which axis is the stacking (out-of-plane) direction: `1`/`2`/`3`,
#'   `"x"`/`"y"`/`"z"`, or (assuming right-anterior-superior orientation)
#'   `"sagittal"`/`"coronal"`/`"axial"`/`"horizontal"`.
#' @return A new `SlicePackage`.
#' @rdname SlicePackage
#' @name SlicePackage_from_image_axis
SlicePackage$from_image_axis <- function(image, axis) {
  if (image$GetDimension() != 3) {
    stop("`image` must be a 3D image.", call. = FALSE)
  }
  axis_index <- .resolve_axis_index(axis)

  size_full <- image$GetSize()
  spacing_full <- image$GetSpacing()
  base <- SliceGeometry$from_image_axis(image, axis_index, coordinate = image$GetOrigin()[axis_index])

  SlicePackage$new(base_slice = base, spacing_k = spacing_full[axis_index], size_k = size_full[axis_index])
}

# Coerce a single item to a SlicePackage: pass a SlicePackage through
# unchanged; wrap a bare SliceGeometry as a single-slice (size_k = 1)
# package. Used so SlicePackageSet can accept either interchangeably.
.as_slice_package <- function(x) {
  if (inherits(x, "SlicePackage")) {
    return(x)
  }
  if (inherits(x, "SliceGeometry")) {
    return(SlicePackage$new(base_slice = x, spacing_k = 1, size_k = 1))
  }
  stop("must be a `SlicePackage` or `SliceGeometry` object", call. = FALSE)
}

#' A collection of SlicePackages, ready for combined intensity extraction
#'
#' @description
#' `SlicePackageSet` aggregates multiple named [SlicePackage] objects (e.g.
#' several different oblique orientations, or several structures of
#' interest) into the single, final long-format table of points from which
#' intensity data are extracted for plotting via the Grammar of Graphics. A
#' bare [SliceGeometry] may be supplied anywhere a `SlicePackage` is expected
#' (in the constructor, `set_package()`, or `set_packages()`) and is
#' automatically wrapped as a single-slice (`size_k = 1`) package.
#'
#' It also handles images with more than 3 dimensions (time, channel,
#' gradient direction, etc.): for each requested combination of non-spatial
#' indices, the corresponding 3D spatial sub-volume is extracted once (and
#' reused across every package that needs it) via `SimpleITK::Extract()`,
#' before each package resamples it in a single call. Images sampled through
#' this class must therefore have at least 3 dimensions — every
#' [SliceGeometry]/[SlicePackage] is fundamentally a plane (or stack of parallel
#' planes) embedded in 3D physical space, so genuinely 1D/2D source images
#' are outside the scope of this sampling model.
#'
#' @export
SlicePackageSet <- R6::R6Class(
  "SlicePackageSet",
  public = list(

    #' @description Create a new `SlicePackageSet`.
    #' @param packages A named list of [SlicePackage] objects (may be empty; add more later with `set_package()`).
    initialize = function(packages = list()) {
      private$set_packages_impl(packages)
    },

    #' @description The underlying named list of `SlicePackage` objects.
    get_packages = function() private$packages_,

    #' @description Names of the packages in this set.
    get_package_names = function() names(private$packages_),

    #' @description Add (or replace) a package.
    #' @param name Single non-empty string identifying the package.
    #' @param package A `SlicePackage` object, or a bare [SliceGeometry]
    #'   (automatically wrapped as a single-slice package).
    set_package = function(name, package) {
      if (!is.character(name) || length(name) != 1 || !nzchar(name)) {
        stop("`name` must be a single non-empty string.", call. = FALSE)
      }
      package <- tryCatch(
        .as_slice_package(package),
        error = function(e) stop("`package` must be a `SlicePackage` or `SliceGeometry` object.", call. = FALSE)
      )
      private$packages_[[name]] <- package
      invisible(self)
    },

    #' @description Remove a package by name.
    #' @param name Single string.
    remove_package = function(name) {
      private$packages_[[name]] <- NULL
      invisible(self)
    },

    #' @description Replace the entire named list of packages at once.
    #' @param packages A named list of [SlicePackage] objects (may be empty).
    set_packages = function(packages) {
      private$set_packages_impl(packages)
      invisible(self)
    },

    #' @description Rename a package without removing/re-adding it.
    #' @param old_name Single string: the package's current name.
    #' @param new_name Single string: its new name. Must not already be in use.
    rename_package = function(old_name, new_name) {
      if (!old_name %in% names(private$packages_)) {
        stop("No package named `", old_name, "` in this set.", call. = FALSE)
      }
      if (new_name %in% names(private$packages_)) {
        stop("A package named `", new_name, "` already exists.", call. = FALSE)
      }
      idx <- which(names(private$packages_) == old_name)
      names(private$packages_)[idx] <- new_name
      invisible(self)
    },

    #' @description Combined sample points across every package, as a tibble
    #'   with columns `package, i, j, k, x, y, z`.
    get_sample_points = function() {
      pkgs <- private$packages_
      if (length(pkgs) == 0) {
        return(tibble::tibble(
          package = character(), i = integer(), j = integer(), k = integer(),
          x = double(), y = double(), z = double()
        ))
      }
      out <- dplyr::bind_rows(lapply(names(pkgs), function(nm) {
        pts <- pkgs[[nm]]$get_sample_points()
        pts$package <- nm
        pts
      }))
      dplyr::select(out, "package", dplyr::everything())
    },

    #' @description Sample intensities for every package — and, for
    #'   higher-dimensional images, every requested combination of
    #'   non-spatial indices — from `image`, combining everything into one
    #'   long-format tibble.
    #' @param image A `SimpleITK` image with dimension 3, 4, or 5. The first
    #'   3 dimensions are treated as spatial; any further dimensions are
    #'   treated as non-spatial (e.g. time, channel) and addressed via `extra_index`.
    #' @param extra_index A named list of 0-indexed index vectors, one per
    #'   non-spatial dimension of `image`, in dimension order (e.g.
    #'   `list(t = 0:9)` for a 4D image with 10 time points). Every
    #'   combination is sampled. Must be an empty list if `image` is 3D.
    #' @param interpolator SimpleITK interpolator name, passed to each package's `sample_intensity()`.
    sample_intensity = function(image, extra_index = list(), interpolator = "sitkLinear") {
      dim <- image$GetDimension()
      if (dim < 3 || dim > 5) {
        stop("`image` must have dimension 3, 4, or 5.", call. = FALSE)
      }
      n_extra <- dim - 3
      if (length(extra_index) != n_extra) {
        stop(sprintf(
          "`image` has %d non-spatial dimension(s); `extra_index` must have exactly %d named entries.",
          n_extra, n_extra
        ), call. = FALSE)
      }

      pkgs <- private$packages_
      if (length(pkgs) == 0) stop("No packages in this SlicePackageSet.", call. = FALSE)

      if (n_extra == 0) {
        combos <- NULL
        n_combo_rows <- 1L
      } else {
        combos <- expand.grid(extra_index, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        n_combo_rows <- nrow(combos)
      }

      subvolume_cache <- new.env(parent = emptyenv())
      results <- vector("list", n_combo_rows * length(pkgs))
      result_i <- 1L

      for (r in seq_len(n_combo_rows)) {
        if (is.null(combos)) {
          sub_image <- image
        } else {
          idx_values <- vapply(names(extra_index), function(nm) as.integer(combos[[nm]][r]), integer(1))
          cache_key <- paste(idx_values, collapse = "_")
          if (!exists(cache_key, envir = subvolume_cache, inherits = FALSE)) {
            extract_size <- as.integer(c(image$GetSize()[1:3], rep(0L, n_extra)))
            extract_index <- as.integer(c(0L, 0L, 0L, idx_values))
            assign(
              cache_key,
              SimpleITK::Extract(image, extract_size, extract_index),
              envir = subvolume_cache
            )
          }
          sub_image <- get(cache_key, envir = subvolume_cache, inherits = FALSE)
        }

        for (nm in names(pkgs)) {
          pts <- pkgs[[nm]]$sample_intensity(sub_image, interpolator = interpolator)
          pts$package <- nm
          if (!is.null(combos)) {
            for (extra_nm in names(extra_index)) pts[[extra_nm]] <- combos[[extra_nm]][r]
          }
          results[[result_i]] <- pts
          result_i <- result_i + 1L
        }
      }

      out <- dplyr::bind_rows(results)
      dplyr::select(out, "package", dplyr::any_of(names(extra_index)), dplyr::everything())
    },

    #' @description Print a short summary of every package in the set.
    #' @param ... Unused; present for compatibility with the generic `print()`.
    print = function(...) {
      pkgs <- private$packages_
      cat("<SlicePackageSet>\n")
      cat("  ", length(pkgs), " package(s)\n", sep = "")
      if (length(pkgs) > 0) {
        name_width <- max(nchar(names(pkgs)))
        for (nm in names(pkgs)) {
          pkg <- pkgs[[nm]]
          cat(
            "    ", formatC(nm, width = -name_width), "  size = ",
            paste(pkg$get_size(), collapse = ", "),
            "  spacing = ", paste(pkg$get_spacing(), collapse = ", "),
            "\n",
            sep = ""
          )
        }
      }
      invisible(self)
    }
  ),
  private = list(
    packages_ = NULL,

    set_packages_impl = function(packages) {
      if (!is.list(packages)) {
        stop("`packages` must be a list of `SlicePackage`/`SliceGeometry` objects.", call. = FALSE)
      }
      if (length(packages) > 0) {
        nms <- names(packages)
        if (is.null(nms) || any(!nzchar(nms))) {
          stop("`packages` must be a fully named list (every element needs a name).", call. = FALSE)
        }
        packages <- tryCatch(
          lapply(packages, .as_slice_package),
          error = function(e) {
            stop("Every element of `packages` must be a `SlicePackage` or `SliceGeometry` object.", call. = FALSE)
          }
        )
      }
      private$packages_ <- packages
    }
  )
)

#' @description
#' Construct a `SlicePackageSet` containing the classic sagittal, coronal,
#' and axial single-slice planes through a reference image, at the given
#' coordinates.
#'
#' @param image A 3D `SimpleITK` image.
#' @param coordinates A fully named list or vector of up to 3 world
#'   coordinates, one per plane to include. Names may be any of the axis
#'   synonyms accepted throughout this package (`1`/`2`/`3`, `"x"`/`"y"`/`"z"`,
#'   `"sagittal"`/`"coronal"`/`"axial"`/`"horizontal"`); any subset of the 3
#'   planes may be supplied, but each plane may only be specified once.
#'   Resulting packages are always named `"sagittal"`, `"coronal"`, `"axial"`,
#'   regardless of which synonym was used to specify them.
#' @return A new `SlicePackageSet`.
#' @rdname SlicePackageSet
#' @name SlicePackageSet_from_orthogonal_triplet
SlicePackageSet$from_orthogonal_triplet <- function(image, coordinates) {
  if (is.null(names(coordinates)) || any(!nzchar(names(coordinates)))) {
    stop("`coordinates` must be a fully named list or vector.", call. = FALSE)
  }
  canonical_names <- c("sagittal", "coronal", "axial")
  seen <- character(0)
  packages <- list()
  for (nm in names(coordinates)) {
    axis_index <- .resolve_axis_index(nm)
    canonical <- canonical_names[axis_index]
    if (canonical %in% seen) {
      stop(
        "`coordinates` specifies the ", canonical, " plane more than once ",
        "(via different synonyms).",
        call. = FALSE
      )
    }
    seen <- c(seen, canonical)
    packages[[canonical]] <- SliceGeometry$from_image_axis(image, axis_index, coordinates[[nm]])
  }
  SlicePackageSet$new(packages)
}
