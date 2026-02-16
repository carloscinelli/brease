# Test prior simulation functions

test_that("sim_brease() returns correct structure", {
  samps <- sim_brease(100)

  expect_true(is.data.frame(samps))
  expect_equal(nrow(samps), 100)
  expect_true(all(c("theta0", "eta_e", "eta_s", "theta1") %in% names(samps)))
})

test_that("sim_brease() theta0 marginal is Beta-distributed", {
  set.seed(42)
  samps <- sim_brease(50000, mu0 = 0.3, n0 = 10)

  # Mean should be close to mu0

  expect_equal(mean(samps$theta0), 0.3, tolerance = 0.01)

  # Variance should be close to Beta variance
  expected_var <- 0.3 * 0.7 / (10 + 1)
  expect_equal(var(samps$theta0), expected_var, tolerance = 0.01)
})

test_that("sim_brease() structural equation holds", {
  samps <- sim_brease(1000)

  expected <- samps$theta0 * (1 - samps$eta_e) + (1 - samps$theta0) * samps$eta_s
  expect_equal(samps$theta1, expected, tolerance = 1e-10)
})

test_that("sim_ib() returns correct structure", {
  samps <- sim_ib(100)

  expect_true(is.data.frame(samps))
  expect_equal(nrow(samps), 100)
  expect_true(all(c("theta0", "theta1") %in% names(samps)))
})

test_that("sim_lt() returns correct structure", {
  samps <- sim_lt(100)

  expect_true(is.data.frame(samps))
  expect_equal(nrow(samps), 100)
  expect_true(all(c("psi", "beta", "theta0", "theta1") %in% names(samps)))
  expect_true(all(samps$theta0 >= 0 & samps$theta0 <= 1))
  expect_true(all(samps$theta1 >= 0 & samps$theta1 <= 1))
})

test_that("log_prior_brease() returns correct value for uniform prior", {
  # Beta(1,1) = Uniform, log density at 0.5 = 0
  lp <- log_prior_brease(0.5, 0.5, 0.5, mu0 = 0.5, n0 = 2,
                          mue = 0.5, ne = 2, mus = 0.5, ns = 2)
  expect_equal(lp, 0, tolerance = 1e-10)
})
