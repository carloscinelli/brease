#' Simulate from the BREASE Prior
#'
#' @description Generates prior predictive samples from the BREASE model.
#'
#' @param n Number of samples to draw.
#' @param mu0 Prior mean for baseline risk (default 0.5).
#' @param n0 Prior sample size for baseline risk (default 2).
#' @param mue Prior mean for efficacy (default 0.5).
#' @param ne Prior sample size for efficacy (default \code{mu0 * n0}).
#' @param mus Prior mean for side-effect risk (default 0.5).
#' @param ns Prior sample size for side-effect risk (default \code{(1 - mu0) * n0}).
#'
#' @return A data.frame with columns: \code{theta0}, \code{eta_e},
#'   \code{eta_s}, \code{theta1}.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Default prior
#' prior_samps <- sim_brease(10000)
#' plot(prior_samps$theta0, prior_samps$theta1,
#'      xlab = expression(theta[0]), ylab = expression(theta[1]),
#'      pch = ".", col = "gray")
#'
#' # Informative prior centered at low risk and moderate efficacy
#' prior_samps <- sim_brease(10000, mu0 = 0.01, n0 = 10,
#'                            mue = 0.5, ne = 5, mus = 0.01, ns = 5)
#'
#' @export
sim_brease <- function(n = 10000,
                       mu0 = 0.5, n0 = 2,
                       mue = 0.5, ne = mu0 * n0,
                       mus = 0.5, ns = (1 - mu0) * n0) {
  out <- list()
  out$theta0 <- rbeta(n, mu0 * n0, (1 - mu0) * n0)
  out$eta_e  <- rbeta(n, mue * ne, (1 - mue) * ne)
  out$eta_s  <- rbeta(n, mus * ns, (1 - mus) * ns)
  out$theta1 <- with(out, theta0 * (1 - eta_e) + (1 - theta0) * eta_s)
  as.data.frame(out)
}


#' Simulate from the Independent Beta Prior
#'
#' @description Generates prior predictive samples from the independent
#'   Beta prior on \eqn{(\theta_0, \theta_1)}.
#'
#' @param n Number of samples to draw.
#' @param mu0 Prior mean for control group risk (default 0.5).
#' @param n0 Prior sample size for control group risk (default 2).
#' @param mu1 Prior mean for treatment group risk (default \code{mu0}).
#' @param n1 Prior sample size for treatment group risk (default \code{n0}).
#'
#' @return A data.frame with columns: \code{theta0}, \code{theta1}.
#'
#' @examples
#' prior_samps <- sim_ib(10000)
#' plot(prior_samps$theta0, prior_samps$theta1, pch = ".", col = "gray")
#'
#' @export
sim_ib <- function(n = 10000,
                   mu0 = 0.5, n0 = 2,
                   mu1 = mu0, n1 = n0) {
  out <- list()
  out$theta0 <- rbeta(n, mu0 * n0, (1 - mu0) * n0)
  out$theta1 <- rbeta(n, mu1 * n1, (1 - mu1) * n1)
  as.data.frame(out)
}


#' Simulate from the Logit-Transform Prior
#'
#' @description Generates prior predictive samples from the logit-transform
#'   (LT) prior, which places normal priors on the log-odds baseline
#'   (\eqn{\beta}) and log-odds ratio (\eqn{\psi}).
#'
#' @param n Number of samples to draw.
#' @param mu_psi Prior mean for log-odds ratio (default 0).
#' @param var_psi Prior variance for log-odds ratio (default 1).
#' @param mu_beta Prior mean for log-odds baseline (default 0).
#' @param var_beta Prior variance for log-odds baseline (default 1).
#'
#' @return A data.frame with columns: \code{psi}, \code{beta},
#'   \code{theta0}, \code{theta1}.
#'
#' @examples
#' prior_samps <- sim_lt(10000)
#' plot(prior_samps$theta0, prior_samps$theta1, pch = ".", col = "gray")
#'
#' @export
sim_lt <- function(n = 10000,
                   mu_psi = 0, var_psi = 1,
                   mu_beta = 0, var_beta = 1) {
  out <- list()
  out$psi    <- rnorm(n, mean = mu_psi, sd = sqrt(var_psi))
  out$beta   <- rnorm(n, mean = mu_beta, sd = sqrt(var_beta))
  out$theta0 <- with(out, exp(beta - psi / 2) / (1 + exp(beta - psi / 2)))
  out$theta1 <- with(out, exp(beta + psi / 2) / (1 + exp(beta + psi / 2)))
  as.data.frame(out)
}


#' Log Prior Density for the BREASE Model
#'
#' @description Evaluates the log prior density of the BREASE model at
#'   given parameter values.
#'
#' @param theta0 Baseline risk.
#' @param eta_e Efficacy.
#' @param eta_s Risk of side effects.
#' @param mu0 Prior mean for baseline risk (default 0.5).
#' @param n0 Prior sample size for baseline risk (default 2).
#' @param mue Prior mean for efficacy (default 0.5).
#' @param ne Prior sample size for efficacy (default 2).
#' @param mus Prior mean for side-effect risk (default 0.5).
#' @param ns Prior sample size for side-effect risk (default 2).
#'
#' @return The log prior density (a scalar).
#'
#' @examples
#' log_prior_brease(0.5, 0.3, 0.1)
#'
#' @export
log_prior_brease <- function(theta0, eta_e, eta_s,
                             mu0 = 0.5, n0 = 2,
                             mue = 0.5, ne = 2,
                             mus = 0.5, ns = 2) {
  dbeta(theta0, mu0 * n0, (1 - mu0) * n0, log = TRUE) +
    dbeta(eta_e, mue * ne, (1 - mue) * ne, log = TRUE) +
    dbeta(eta_s, mus * ns, (1 - mus) * ns, log = TRUE)
}
