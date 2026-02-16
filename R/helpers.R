# Internal helper functions for the brease package
# These are not exported.

#' Compute treatment effect statistics from posterior samples
#'
#' Computes risk ratio (RR), efficacy fraction (EF), odds ratio (OR),
#' and risk difference (RD) from posterior samples of theta0 and theta1.
#'
#' @param samps A list or data.frame containing \code{theta0} and \code{theta1}.
#' @return The input with additional columns: \code{rr}, \code{ef}, \code{or}, \code{rd}.
#' @keywords internal
compute_stats <- function(samps) {
  samps <- samps[names(samps) %in% c("theta0", "theta1", "eta_e", "eta_s")]
  samps$rr <- with(samps, theta1 / theta0)
  samps$ef <- with(samps, 1 - theta1 / theta0)
  samps$or <- with(samps, (theta1 / (1 - theta1)) / (theta0 / (1 - theta0)))
  samps$rd <- with(samps, theta1 - theta0)
  return(samps)
}

#' Compute posterior summaries
#'
#' Computes mean, median, and credible interval for each column in a
#' data.frame of posterior samples.
#'
#' @param samps A data.frame of posterior samples.
#' @param med Quantile for the posterior median (default 0.5).
#' @param lw Lower quantile for the credible interval.
#' @param up Upper quantile for the credible interval.
#' @return A matrix with rows for each variable and columns: mean, med, lw, up.
#' @keywords internal
compute_summaries <- function(samps, med = 0.5, lw = 0.025, up = 0.975) {
  f <- function(x) c(mean(x), quantile(x, c(med, lw, up)))
  post_summaries <- t(sapply(samps, f))
  colnames(post_summaries) <- c("mean", "med", "lw", "up")
  return(post_summaries)
}

#' Structural equation for theta1
#'
#' Computes the treatment group risk from baseline risk, efficacy,
#' and side-effect risk.
#'
#' @param theta0 Baseline risk.
#' @param eta_e Efficacy.
#' @param eta_s Risk of side effects.
#' @return Treatment group risk: \code{theta0 * (1 - eta_e) + (1 - theta0) * eta_s}.
#' @keywords internal
f_theta1 <- function(theta0, eta_e, eta_s) {
  theta0 * (1 - eta_e) + (1 - theta0) * eta_s
}


# IB helpers (internal) ---------------------------------------------------

#' Independent Beta log marginal likelihood under H0
#' @keywords internal
lml_ib_h0 <- function(y0, y1, N0, N1, a0 = 1, b0 = 1) {
  lchoose(N0, y0) + lchoose(N1, y1) +
    lbeta(y0 + y1 + a0, N0 + N1 - (y0 + y1) + b0) - lbeta(a0, b0)
}

#' Independent Beta log marginal likelihood under H1
#' @keywords internal
lml_ib_h1 <- function(y0, y1, N0, N1, a0 = 1, b0 = 1, a1 = 1, b1 = 1) {
  lchoose(N0, y0) + lchoose(N1, y1) +
    lbeta(y0 + a0, N0 - y0 + b0) - lbeta(a0, b0) +
    lbeta(y1 + a1, N1 - y1 + b1) - lbeta(a1, b1)
}

#' Independent Beta analytic Bayes factor (log scale)
#' @keywords internal
bf_ib_analytic <- function(y0, y1, N0, N1, a = 1, b = 1) {
  post0 <- lbeta(2 * a - 1, 2 * b - 1) - 2 * lbeta(a, b)
  post1 <- lbeta(2 * a + y0 + y1 - 1, 2 * b + N0 - y0 + N1 - y1 - 1) -
    lbeta(a + y0, b + N0 - y0) - lbeta(a + y1, b + N1 - y1)
  post0 - post1
}

# Input validation helpers ------------------------------------------------

#' Check that a value is a non-negative integer (count)
#' @keywords internal
check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || x < 0 || x != round(x)) {
    stop(sprintf("'%s' must be a non-negative integer.", name), call. = FALSE)
  }
}

#' Check that a value is a positive integer
#' @keywords internal
check_positive_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || x <= 0 || x != round(x)) {
    stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
  }
}

#' Check that a value is a probability (between 0 and 1)
#' @keywords internal
check_prob <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || x < 0 || x > 1) {
    stop(sprintf("'%s' must be a number between 0 and 1.", name), call. = FALSE)
  }
}

#' Check that a value is a positive number
#' @keywords internal
check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || x <= 0) {
    stop(sprintf("'%s' must be a positive number.", name), call. = FALSE)
  }
}

#' Validate BREASE data inputs
#' @keywords internal
check_data <- function(y0, y1, N0, N1) {
  check_count(y0, "y0")
  check_count(y1, "y1")
  check_positive_int(N0, "N0")
  check_positive_int(N1, "N1")
  if (y0 > N0) stop("'y0' cannot be greater than 'N0'.", call. = FALSE)
  if (y1 > N1) stop("'y1' cannot be greater than 'N1'.", call. = FALSE)
}

#' Validate BREASE prior hyperparameters
#' @keywords internal
check_prior <- function(mu0, n0, mue, ne, mus, ns) {
  check_prob(mu0, "mu0")
  check_positive(n0, "n0")
  check_prob(mue, "mue")
  check_positive(ne, "ne")
  check_prob(mus, "mus")
  check_positive(ns, "ns")
}
