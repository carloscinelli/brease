# Test marginal likelihood and Bayes factor functions

test_that("lml_brease() returns a finite number", {
  lml <- lml_brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
  expect_true(is.numeric(lml))
  expect_true(is.finite(lml))
})

test_that("lml_brease_h0() returns a finite number", {
  lml <- lml_brease_h0(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
  expect_true(is.numeric(lml))
  expect_true(is.finite(lml))
})

test_that("lml_brease_mono() returns a finite number", {
  lml <- lml_brease_mono(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)
  expect_true(is.numeric(lml))
  expect_true(is.finite(lml))
})

test_that("bf_brease() returns correct log Bayes factor", {
  lbf01 <- bf_brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                     mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
                     mus = 0.3, ns = 1)

  # Should equal lml_h0 - lml_h1
  lml_h0 <- lml_brease_h0(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
  lml_h1 <- lml_brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                        mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
                        mus = 0.3, ns = 1)
  expect_equal(lbf01, lml_h0 - lml_h1, tolerance = 1e-10)
})

test_that("bf_brease_mono() is consistent", {
  lbf01 <- bf_brease_mono(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965,
                          mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1)

  lml_h0 <- lml_brease_h0(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)
  lml_h1 <- lml_brease_mono(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965,
                             mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1)
  expect_equal(lbf01, lml_h0 - lml_h1, tolerance = 1e-10)
})

test_that("H0 marginal likelihood is always less than H1 when data shows clear effect", {
  # COVID data shows very strong effect
  lml_h0 <- lml_brease_h0(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)
  lml_h1 <- lml_brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)

  # H1 should fit the data better
  expect_true(lml_h1 > lml_h0)
})

test_that("lml_brease_h0() with uniform prior matches IB formula", {
  y0 <- 26; y1 <- 10; N0 <- 11034; N1 <- 11037

  # BREASE H0 with mu0 = 0.5, n0 = 2 gives Beta(1,1) = Uniform
  brease_h0 <- lml_brease_h0(y0, y1, N0, N1, mu0 = 0.5, n0 = 2)

  # IB H0 with a = 1, b = 1
  ib_h0 <- lchoose(N0, y0) + lchoose(N1, y1) +
    lbeta(y0 + y1 + 1, N0 + N1 - (y0 + y1) + 1) - lbeta(1, 1)

  expect_equal(brease_h0, ib_h0, tolerance = 1e-10)
})
