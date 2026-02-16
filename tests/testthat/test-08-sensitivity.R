# Test sensitivity analysis functions

test_that("bf_seq() returns correct structure", {
  sens <- bf_seq(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                 mus_seq = seq(0.1, 0.9, length = 5),
                 mue_seq = seq(0.1, 0.9, length = 5))

  expect_true(is.list(sens))
  expect_true("mus_seq" %in% names(sens))
  expect_true("mue_seq" %in% names(sens))
  expect_true("log_bf_seq" %in% names(sens))
  expect_equal(dim(sens$log_bf_seq), c(5, 5))
})

test_that("bf_seq() log BFs are finite", {
  sens <- bf_seq(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                 mus_seq = seq(0.1, 0.9, length = 3),
                 mue_seq = seq(0.1, 0.9, length = 3))

  expect_true(all(is.finite(sens$log_bf_seq)))
})
