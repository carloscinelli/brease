# Test posterior sampling functions

test_that("sample_brease() returns correct structure", {
  samps <- sample_brease(100, y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  expect_true(is.list(samps))
  expect_equal(length(samps$theta0), 100)
  expect_equal(length(samps$theta1), 100)
  expect_equal(length(samps$eta_e), 100)
  expect_equal(length(samps$eta_s), 100)
})

test_that("sample_brease() samples are in valid range", {
  samps <- sample_brease(1000, y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  expect_true(all(samps$theta0 >= 0 & samps$theta0 <= 1))
  expect_true(all(samps$theta1 >= 0 & samps$theta1 <= 1))
  expect_true(all(samps$eta_e >= 0 & samps$eta_e <= 1))
  expect_true(all(samps$eta_s >= 0 & samps$eta_s <= 1))
})

test_that("sample_brease() satisfies structural equation", {
  samps <- sample_brease(1000, y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  expected_theta1 <- samps$theta0 * (1 - samps$eta_e) +
    (1 - samps$theta0) * samps$eta_s
  expect_equal(samps$theta1, expected_theta1, tolerance = 1e-10)
})

test_that("sample_brease_mono() has zero side effects", {
  samps <- sample_brease_mono(100, y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)

  expect_true(all(samps$eta_s == 0))
  expect_equal(samps$theta1, samps$theta0 * (1 - samps$eta_e), tolerance = 1e-10)
})

test_that("gibbs_brease() returns data.frame with correct columns", {
  samps <- gibbs_brease(200, y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                         burn_in = 100)

  expect_true(is.data.frame(samps))
  expect_equal(nrow(samps), 200)
  expect_true(all(c("theta0", "theta1", "eta_e", "eta_s") %in% names(samps)))
})

test_that("exact sampler and Gibbs sampler produce similar posteriors", {
  set.seed(42)
  n <- 5000
  y0 <- 26; y1 <- 10; N0 <- 11034; N1 <- 11037

  exact <- sample_brease(n, y0, y1, N0, N1,
                          mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
                          mus = 0.3, ns = 1)
  gibbs <- gibbs_brease(n, y0, y1, N0, N1,
                         mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
                         mus = 0.3, ns = 1, burn_in = 2000)

  # Compare posterior means (should be within Monte Carlo error)
  expect_equal(mean(exact$theta0), mean(gibbs$theta0), tolerance = 0.01)
  expect_equal(mean(exact$theta1), mean(gibbs$theta1), tolerance = 0.01)
})
