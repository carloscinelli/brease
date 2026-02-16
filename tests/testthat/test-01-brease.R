# Test the main brease() function

test_that("brease() returns correct S3 structure", {
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  expect_s3_class(result, "brease")
  expect_named(result, c("post.samples", "post.summaries", "bayes.factors",
                          "call", "data", "prior", "info"))
  expect_true(is.data.frame(result$post.samples))
  expect_true(is.matrix(result$post.summaries))
  expect_true(is.list(result$bayes.factors))
  expect_true("bf10" %in% names(result$bayes.factors))
})

test_that("brease() posterior samples have correct columns", {
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037, n_samples = 100)

  expect_true("theta0" %in% names(result$post.samples))
  expect_true("theta1" %in% names(result$post.samples))
  expect_true("rr"     %in% names(result$post.samples))
  expect_true("ef"     %in% names(result$post.samples))
  expect_true("or"     %in% names(result$post.samples))
  expect_true("rd"     %in% names(result$post.samples))
  expect_equal(nrow(result$post.samples), 100)
})

test_that("brease() posterior summaries have correct format", {
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037, n_samples = 100)

  ps <- result$post.summaries
  expect_equal(colnames(ps), c("mean", "med", "lw", "up"))
  expect_true("theta0" %in% rownames(ps))
  expect_true("rr"     %in% rownames(ps))
})

test_that("brease() with mono = TRUE works", {
  result <- brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965, mono = TRUE,
                   n_samples = 100)

  expect_s3_class(result, "brease")
  expect_true(result$info$mono)
  expect_true(all(result$post.samples$eta_s == 0))
})

test_that("brease() Gibbs sampler works", {
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                   sampler = "gibbs", n_samples = 100, burn_in = 100)

  expect_s3_class(result, "brease")
  expect_equal(result$info$sampler, "gibbs")
})

test_that("print.brease() works without error", {
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037, n_samples = 100)
  expect_output(print(result), "Bayesian Analysis")
})

test_that("summary.brease() works without error", {
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037, n_samples = 100)
  expect_output(summary(result), "Posterior Summaries")
})

test_that("brease() validates inputs", {
  expect_error(brease(y0 = -1, y1 = 10, N0 = 100, N1 = 100))
  expect_error(brease(y0 = 10, y1 = 10, N0 = 5, N1 = 100))
  expect_error(brease(y0 = 10, y1 = 10, N0 = 100, N1 = 100, mu0 = 2))
})
