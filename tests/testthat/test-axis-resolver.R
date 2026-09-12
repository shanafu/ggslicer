test_that(".resolve_axis_index() accepts numeric, cartesian, and anatomical forms", {
  expect_equal(.resolve_axis_index(1), 1L)
  expect_equal(.resolve_axis_index("1"), 1L)
  expect_equal(.resolve_axis_index("x"), 1L)
  expect_equal(.resolve_axis_index("X"), 1L)
  expect_equal(.resolve_axis_index("sagittal"), 1L)
  expect_equal(.resolve_axis_index("Sagittal"), 1L)

  expect_equal(.resolve_axis_index(2), 2L)
  expect_equal(.resolve_axis_index("y"), 2L)
  expect_equal(.resolve_axis_index("coronal"), 2L)

  expect_equal(.resolve_axis_index(3), 3L)
  expect_equal(.resolve_axis_index("z"), 3L)
  expect_equal(.resolve_axis_index("axial"), 3L)
  expect_equal(.resolve_axis_index("horizontal"), 3L)
  expect_equal(.resolve_axis_index("Horizontal"), 3L)
})

test_that(".resolve_axis_index() rejects invalid input", {
  expect_error(.resolve_axis_index(0), "1, 2, or 3")
  expect_error(.resolve_axis_index(4), "1, 2, or 3")
  expect_error(.resolve_axis_index(c(1, 2)), "1, 2, or 3")
  expect_error(.resolve_axis_index("bogus"), "Invalid `axis`")
  expect_error(.resolve_axis_index(TRUE), "must be a single value")
})
