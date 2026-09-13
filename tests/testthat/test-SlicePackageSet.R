make_package <- function(direction_i = c(1, 0, 0), direction_j = c(0, 1, 0), size = c(3, 3), spacing_k = 1, size_k = 2) {
  base <- SliceGeometry$new(
    origin = c(0, 0, 0), direction_i = direction_i, direction_j = direction_j,
    spacing = c(1, 1), size = size
  )
  SlicePackage$new(base_slice = base, spacing_k = spacing_k, size_k = size_k)
}

test_that("constructor validates its inputs", {
  expect_error(SlicePackageSet$new(list(a = 1, b = 2)), "SlicePackage")
  expect_error(SlicePackageSet$new(list(a = make_package(), 2)), "fully named")
  expect_error(SlicePackageSet$new("not a list"), "list")
  expect_silent(SlicePackageSet$new())
  expect_silent(SlicePackageSet$new(list(a = make_package())))
})

test_that("set_package()/remove_package()/get_package_names() manage the collection", {
  sset <- SlicePackageSet$new()
  expect_equal(sset$get_package_names(), NULL)

  sset$set_package("axial", make_package())
  expect_equal(sset$get_package_names(), "axial")

  expect_error(sset$set_package("axial", "not a package"), "SlicePackage")
  expect_error(sset$set_package(1, make_package()), "non-empty string")

  sset$set_package("sagittal", make_package(direction_i = c(0, 1, 0), direction_j = c(0, 0, 1)))
  expect_setequal(sset$get_package_names(), c("axial", "sagittal"))

  sset$remove_package("axial")
  expect_equal(sset$get_package_names(), "sagittal")
})

test_that("get_sample_points() combines packages with a package column", {
  sset <- SlicePackageSet$new(list(
    axial = make_package(size_k = 4),
    sagittal = make_package(direction_i = c(0, 1, 0), direction_j = c(0, 0, 1), size_k = 2)
  ))
  pts <- sset$get_sample_points()

  expect_equal(names(pts)[1], "package")
  expect_setequal(unique(pts$package), c("axial", "sagittal"))
  expect_equal(nrow(pts), 3 * 3 * 4 + 3 * 3 * 2)
})

test_that("get_sample_points() on an empty set returns a typed, 0-row tibble", {
  sset <- SlicePackageSet$new()
  pts <- sset$get_sample_points()

  expect_equal(nrow(pts), 0)
  expect_named(pts, c("package", "i", "j", "k", "x", "y", "z"))
})

test_that("sample_intensity() works directly on a 3D image with no extra_index", {
  skip_if_not_installed("SimpleITK")

  sset <- SlicePackageSet$new(list(only = make_package(size = c(5, 5), size_k = 10)))
  src <- filled_sitk_image(c(5, 5, 10), function(idx) 100 * idx[3] + 10 * idx[1] + idx[2])

  out <- sset$sample_intensity(src)
  expect_equal(nrow(out), 5 * 5 * 10)
  row <- out[out$i == 1 & out$j == 2 & out$k == 3, ]
  expect_equal(row$intensity, 100 * 3 + 10 * 1 + 2)
})

test_that("sample_intensity() handles a 4D image via extract-and-cache, matching direct pixel lookups", {
  skip_if_not_installed("SimpleITK")

  img4d <- filled_sitk_image(
    c(5, 5, 10, 3),
    function(idx) 1000 * idx[4] + 100 * idx[3] + 10 * idx[1] + idx[2]
  )

  sset <- SlicePackageSet$new(list(
    axial = make_package(size = c(5, 5), size_k = 10),
    sagittal = make_package(direction_i = c(0, 1, 0), direction_j = c(0, 0, 1), size = c(5, 10), size_k = 5)
  ))

  out <- sset$sample_intensity(img4d, extra_index = list(t = 0:2))

  expect_equal(nrow(out), (5 * 5 * 10 + 5 * 10 * 5) * 3)
  expect_true(all(c("package", "t") %in% names(out)))

  chk <- out[out$package == "axial" & out$t == 2 & out$i == 2 & out$j == 1 & out$k == 3, ]
  expect_equal(nrow(chk), 1)
  expect_equal(chk$intensity, img4d$GetPixel(c(2L, 1L, 3L, 2L)))
})

test_that("sample_intensity() handles a 5D image (two non-spatial dimensions)", {
  skip_if_not_installed("SimpleITK")

  img5d <- filled_sitk_image(
    c(2, 2, 2, 2, 3),
    function(idx) 1000 * idx[5] + 100 * idx[4] + 10 * idx[1] + idx[2] + idx[3]
  )
  sset <- SlicePackageSet$new(list(only = make_package(size = c(2, 2), size_k = 2)))

  out <- sset$sample_intensity(img5d, extra_index = list(t = 0:1, c = 0:2))

  expect_equal(nrow(out), 2 * 2 * 2 * 2 * 3)
  chk <- out[out$t == 1 & out$c == 2 & out$i == 1 & out$j == 0 & out$k == 1, ]
  expect_equal(nrow(chk), 1)
  expect_equal(chk$intensity, img5d$GetPixel(c(1L, 0L, 1L, 1L, 2L)))
})

test_that("sample_intensity() validates image dimension against extra_index", {
  skip_if_not_installed("SimpleITK")

  sset <- SlicePackageSet$new(list(only = make_package()))
  src3d <- filled_sitk_image(c(3, 3, 2), function(idx) 0)
  img4d <- filled_sitk_image(c(2, 2, 2, 2), function(idx) 0)

  expect_error(sset$sample_intensity(src3d, extra_index = list(t = 0:1)), "0 non-spatial")
  expect_error(sset$sample_intensity(img4d), "1 non-spatial")
  expect_error(sset$sample_intensity(img4d, extra_index = list(t = 0:1, c = 0:1)), "1 non-spatial")
})

test_that("sample_intensity() errors when the set has no packages", {
  skip_if_not_installed("SimpleITK")

  sset <- SlicePackageSet$new()
  src3d <- filled_sitk_image(c(3, 3, 2), function(idx) 0)
  expect_error(sset$sample_intensity(src3d), "No packages")
})

test_that("set_packages() bulk-replaces the whole collection", {
  sset <- SlicePackageSet$new(list(a = make_package()))
  sset$set_packages(list(x = make_package(), y = make_package(size_k = 3)))

  expect_setequal(sset$get_package_names(), c("x", "y"))

  sset$set_packages(list())
  expect_null(sset$get_package_names())

  expect_error(sset$set_packages(list(1, 2)), "fully named")
})

test_that("rename_package() renames without disturbing the underlying package", {
  pkg <- make_package()
  sset <- SlicePackageSet$new(list(old = pkg))
  sset$rename_package("old", "new")

  expect_equal(sset$get_package_names(), "new")
  expect_identical(sset$get_packages()[["new"]], pkg)

  expect_error(sset$rename_package("missing", "z"), "No package named")

  sset$set_package("another", make_package())
  expect_error(sset$rename_package("new", "another"), "already exists")
})

test_that("print() reports each package's name, size, and spacing", {
  sset <- SlicePackageSet$new(list(axial = make_package(), sagittal = make_package(size_k = 3)))
  expect_output(print(sset), "<SlicePackageSet>")
  expect_output(print(sset), "axial")
  expect_output(print(sset), "sagittal")
  expect_output(print(sset), "size = ")

  expect_output(print(SlicePackageSet$new()), "0 package")
})

test_that("a bare SliceGeometry is accepted wherever a SlicePackage is expected", {
  slice <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))

  sset <- SlicePackageSet$new(list(single = slice))
  wrapped <- sset$get_packages()[["single"]]
  expect_s3_class(wrapped, "SlicePackage")
  expect_equal(wrapped$get_size_k(), 1L)
  expect_equal(wrapped$get_base_slice()$get_origin(), slice$get_origin())

  sset$set_package("another", slice)
  expect_s3_class(sset$get_packages()[["another"]], "SlicePackage")

  sset$set_packages(list(mixed_a = slice, mixed_b = make_package()))
  expect_true(all(vapply(sset$get_packages(), inherits, logical(1), what = "SlicePackage")))

  expect_error(SlicePackageSet$new(list(bad = "not a slice or package")), "SlicePackage.*SliceGeometry")
})

test_that("SlicePackageSet$from_orthogonal_triplet() builds the classic 3-plane view", {
  skip_if_not_installed("SimpleITK")

  img <- filled_sitk_image(c(10, 10, 10), function(idx) 0, origin = c(-5, -5, -5), spacing = c(1, 1, 1))

  triplet <- SlicePackageSet$from_orthogonal_triplet(img, list(x = 0, coronal = 1, horizontal = -2))
  expect_setequal(triplet$get_package_names(), c("sagittal", "coronal", "axial"))
  expect_equal(triplet$get_packages()[["sagittal"]]$get_size_k(), 1L)

  partial <- SlicePackageSet$from_orthogonal_triplet(img, list(axial = 0))
  expect_equal(partial$get_package_names(), "axial")

  expect_error(
    SlicePackageSet$from_orthogonal_triplet(img, list(x = 0, sagittal = 1)),
    "more than once"
  )
  expect_error(SlicePackageSet$from_orthogonal_triplet(img, list(0, 1)), "fully named")
})

test_that("SlicePackageSet$from_slice_packages() builds a set from an unnamed list or a single package", {
  base <- SliceGeometry$new(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1), c(3, 3))
  pkg1 <- SlicePackage$new(base_slice = base, spacing_k = 1, size_k = 2)
  pkg2 <- SlicePackage$new(base_slice = base, spacing_k = 1, size_k = 3)

  sset <- SlicePackageSet$from_slice_packages(list(pkg1, pkg2))
  expect_equal(sset$get_package_names(), c("package_1", "package_2"))
  expect_equal(nrow(sset$get_sample_points()), 9 * 2 + 9 * 3)

  # a single SlicePackage, not wrapped in a list
  sset_single <- SlicePackageSet$from_slice_packages(pkg1)
  expect_equal(sset_single$get_package_names(), "package_1")
  expect_equal(nrow(sset_single$get_sample_points()), 9 * 2)

  # a single bare SliceGeometry, auto-wrapped as size_k = 1
  sset_slice <- SlicePackageSet$from_slice_packages(base)
  expect_equal(sset_slice$get_packages()[["package_1"]]$get_size_k(), 1L)

  # a mixed list of SlicePackage and bare SliceGeometry
  sset_mixed <- SlicePackageSet$from_slice_packages(list(pkg1, base))
  expect_equal(sset_mixed$get_package_names(), c("package_1", "package_2"))
  expect_equal(sset_mixed$get_packages()[["package_2"]]$get_size_k(), 1L)

  expect_error(SlicePackageSet$from_slice_packages("not a package"), "SlicePackage.*SliceGeometry")
})
