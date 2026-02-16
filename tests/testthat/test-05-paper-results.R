# Test replication of paper results

test_that("Aspirin trial: BF10 shows evidence for treatment effect", {
  # From Table in Irons & Cinelli (2025)
  # With mue = mus = 0.3, ne = ns = 1, the aspirin trial
  # should show moderate evidence for treatment effect
  result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
                   mue = 0.3, mus = 0.3, ne = 1, ns = 1)

  # BF10 should be > 1 (evidence for effect)
  expect_true(result$bayes.factors$bf10 > 1)
})

test_that("COVID trial: strong evidence for vaccine efficacy", {
  result <- brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965,
                   mue = 0.3, mus = 0.3, ne = 1, ns = 1)

  # Very strong evidence for treatment effect
  expect_true(result$bayes.factors$bf10 > 100)

  # Efficacy fraction should be very high (> 90%)
  expect_true(result$post.summaries["ef", "mean"] > 0.90)
})

test_that("COVID trial: monotone model gives strong evidence", {
  result <- brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965,
                   mue = 0.3, ne = 1, mono = TRUE)

  # Monotone model should also show very strong evidence
  expect_true(result$bayes.factors$bf10 > 100)
})

test_that("Aspirin trial: IB comparison produces reasonable results", {
  result <- brease_ib(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)

  # Risk ratio should be < 1 (aspirin reduces MI risk)
  expect_true(result$post.summaries["rr", "mean"] < 1)
})
