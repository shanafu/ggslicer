test_that("quantile method returns exactly n interior quantiles", {
  x <- 1:100
  levels <- suggest_contour_levels(data.frame(value = x), method = "quantile", n = 5)
  expected <- unname(quantile(x, probs = (1:5) / 6))
  expect_equal(as.numeric(levels), as.numeric(expected))
  expect_length(levels, 5)
})

test_that("quantile method respects a custom n", {
  x <- 1:100
  levels <- suggest_contour_levels(data.frame(value = x), method = "quantile", n = 3)
  expect_length(levels, 3)
})

test_that("troughs method finds the real boundary of a well-separated bimodal distribution", {
  for (seed in 1:5) {
    set.seed(seed)
    bimodal <- c(rnorm(1000, mean = 0, sd = 1), rnorm(1000, mean = 20, sd = 1))
    levels <- suggest_contour_levels(data.frame(value = bimodal), method = "troughs", min_n = 1)
    expect_true(all(levels > 3 & levels < 17), info = paste("seed", seed))
  }
})

test_that("troughs method pads a unimodal distribution to exactly min_n levels", {
  for (seed in 1:5) {
    set.seed(seed)
    unimodal <- rnorm(1000)
    levels <- suggest_contour_levels(data.frame(value = unimodal), method = "troughs", min_n = 3)
    expect_length(levels, 3)
  }
})

test_that("troughs method finds multiple real boundaries for a trimodal distribution", {
  set.seed(1)
  trimodal <- c(rnorm(1000, 0, 1), rnorm(1000, 20, 1), rnorm(1000, 40, 1))
  levels <- suggest_contour_levels(data.frame(value = trimodal), method = "troughs", min_n = 1)
  expect_gte(length(levels), 2)
  expect_true(any(levels > 3 & levels < 17))
  expect_true(any(levels > 23 & levels < 37))
})

test_that("both a data frame and a SimpleITK image/path work as input", {
  skip_if_not_installed("SimpleITK")
  set.seed(1)
  img <- filled_sitk_image(c(10, 10, 10), function(idx) sum(idx) + runif(1))

  levels_img <- suggest_contour_levels(img, method = "quantile", n = 3)
  expect_length(levels_img, 3)

  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path))
  SimpleITK::WriteImage(img, path)
  levels_path <- suggest_contour_levels(path, method = "quantile", n = 3)
  expect_equal(levels_img, levels_path)
})

test_that("a column name matching discrete_data_names() errors clearly", {
  df <- data.frame(mask = c(0, 1, 0, 1))
  expect_error(suggest_contour_levels(df, column = "mask"), "discrete")
})

test_that("data-driven discreteness check catches an oddly-named discrete column", {
  df <- data.frame(oddname = rep(0:3, 25))
  expect_error(suggest_contour_levels(df, column = "oddname"), "discrete")
})

test_that("a genuinely continuous, oddly-named column is not flagged as discrete", {
  set.seed(1)
  df <- data.frame(oddname = rnorm(1000))
  expect_no_error(suggest_contour_levels(df, column = "oddname"))
})

test_that("suggest_contour_levels() validates its arguments", {
  df <- data.frame(value = rnorm(100))
  expect_error(suggest_contour_levels(df, method = "bogus"))
  expect_error(suggest_contour_levels(df, n = 0))
  expect_error(suggest_contour_levels(df, n = -1))
  expect_error(suggest_contour_levels(df, method = "troughs", min_n = 0))
  expect_error(suggest_contour_levels(df, column = "missing"), "no column")
})

test_that("suggested levels work directly as slice_contours()'s levels argument", {
  skip_if_not_installed("SimpleITK")
  set.seed(1)
  img <- filled_sitk_image(c(20, 20, 20), function(idx) {
    sqrt((idx[1] - 10)^2 + (idx[2] - 10)^2)
  })
  levels <- suggest_contour_levels(img, method = "quantile", n = 3)
  out <- slice_contours(img, axis = "z", coordinate = 10, levels = levels)
  expect_true(nrow(out) > 0)
})
