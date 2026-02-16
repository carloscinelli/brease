# Test bounds on probabilities of causation

test_that("brease_bounds() returns correct structure", {
  bounds <- brease_bounds(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  expect_true(is.list(bounds))
  expect_named(bounds, c("eta_s_lower", "eta_s_upper", "eta_s_mid",
                          "eta_e_lower", "eta_e_upper", "eta_e_mid"))
})

test_that("brease_bounds() bounds are valid probabilities", {
  bounds <- brease_bounds(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  expect_true(bounds$eta_s_lower >= 0 && bounds$eta_s_lower <= 1)
  expect_true(bounds$eta_s_upper >= 0 && bounds$eta_s_upper <= 1)
  expect_true(bounds$eta_e_lower >= 0 && bounds$eta_e_lower <= 1)
  expect_true(bounds$eta_e_upper >= 0 && bounds$eta_e_upper <= 1)
  expect_true(bounds$eta_s_lower <= bounds$eta_s_upper)
  expect_true(bounds$eta_e_lower <= bounds$eta_e_upper)
})

test_that("brease_bounds() midpoints are averages of bounds", {
  bounds <- brease_bounds(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)

  expect_equal(bounds$eta_s_mid,
               (bounds$eta_s_lower + bounds$eta_s_upper) / 2)
  expect_equal(bounds$eta_e_mid,
               (bounds$eta_e_lower + bounds$eta_e_upper) / 2)
})

test_that("brease_bounds() validates inputs", {
  expect_error(brease_bounds(y0 = -1, y1 = 10, N0 = 100, N1 = 100))
  expect_error(brease_bounds(y0 = 10, y1 = 10, N0 = 5, N1 = 100))
})
