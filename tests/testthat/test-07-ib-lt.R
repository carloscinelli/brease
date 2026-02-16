# Test IB and LT methods

test_that("brease_ib() returns correct structure", {
  result <- brease_ib(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                      n_samples = 100)

  expect_s3_class(result, "brease")
  expect_equal(result$info$method, "ib")
  expect_true("bf10" %in% names(result$bayes.factors))
})

test_that("brease_ib() with custom priors works", {
  result <- brease_ib(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                      a = c(17, 17), b = c(17, 17), n_samples = 100)

  expect_s3_class(result, "brease")
})

test_that("brease_ib() conjugate posterior is correct", {
  set.seed(42)
  y0 <- 26; N0 <- 11034; a0 <- 1; b0 <- 1

  result <- brease_ib(y0 = y0, y1 = 10, N0 = N0, N1 = 11037,
                      n_samples = 50000)

  # Posterior mean for theta0 should be (y0 + a) / (N0 + a + b)
  expected_mean <- (y0 + a0) / (N0 + a0 + b0)
  expect_equal(mean(result$post.samples$theta0), expected_mean, tolerance = 0.001)
})

test_that("brease_lt() fails gracefully without abtest", {
  # Skip this test if abtest is available (we're testing the error)
  skip_if(requireNamespace("abtest", quietly = TRUE),
          "abtest is installed; cannot test missing package error")
  expect_error(brease_lt(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037),
               "abtest")
})

test_that("print works for IB results", {
  result <- brease_ib(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                      n_samples = 100)
  expect_output(print(result), "Independent Beta")
})
