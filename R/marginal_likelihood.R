#' Log Marginal Likelihood for the BREASE Model (H1)
#'
#' @description Computes the log marginal likelihood under the alternative
#'   hypothesis (treatment has an effect) for the BREASE model. This is an
#'   exact analytic calculation.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#' @param mu0 Prior mean for baseline risk (default 0.5).
#' @param n0 Prior sample size for baseline risk (default 2).
#' @param mue Prior mean for efficacy (default 0.5).
#' @param ne Prior sample size for efficacy (default 2).
#' @param mus Prior mean for side-effect risk (default 0.5).
#' @param ns Prior sample size for side-effect risk (default 2).
#'
#' @return The log marginal likelihood (a scalar).
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial
#' lml_brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
#'            mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1, mus = 0.3, ns = 1)
#'
#' @export
lml_brease <- function(y0, y1, N0, N1,
                       mu0 = 0.5, n0 = 2,
                       mue = 0.5, ne = 2,
                       mus = 0.5, ns = 2) {
  N <- N0 + N1

  lml <- lchoose(N0, y0) + lchoose(N1, y1) -
    lbeta(mu0 * n0, (1 - mu0) * n0) -
    lbeta(mue * ne, (1 - mue) * ne) -
    lbeta(mus * ns, (1 - mus) * ns)

  J <- (0:y1) %o% rep(1, N1 - y1 + 1)
  K <- rep(1, y1 + 1) %o% (0:(N1 - y1))

  summands <- lchoose(y1, J) + lchoose(N1 - y1, K) +
    lbeta(y0 + J + K + mu0 * n0, N - (y0 + J + K) + (1 - mu0) * n0) +
    lbeta(K + mue * ne, J + (1 - mue) * ne) +
    lbeta(y1 - J + mus * ns, N1 - y1 - K + (1 - mus) * ns)

  max_log <- max(summands)
  lml <- lml + log(sum(exp(summands - max_log))) + max_log
  return(lml)
}


#' Log Marginal Likelihood under the Sharp Null (H0)
#'
#' @description Computes the log marginal likelihood under the sharp null
#'   hypothesis \eqn{\theta_0 = \theta_1} for the BREASE model.
#'
#' @inheritParams lml_brease
#'
#' @return The log marginal likelihood under H0 (a scalar).
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' lml_brease_h0(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#'
#' @export
lml_brease_h0 <- function(y0, y1, N0, N1, mu0 = 0.5, n0 = 2) {
  a0 <- mu0 * n0
  b0 <- (1 - mu0) * n0
  lml_ib_h0(y0, y1, N0, N1, a0, b0)
}


#' Log Marginal Likelihood for the Monotone BREASE Model (H1)
#'
#' @description Computes the log marginal likelihood under the alternative
#'   hypothesis for the monotone (no-harm) BREASE model, where
#'   \eqn{\eta_s = 0}.
#'
#' @inheritParams lml_brease
#'
#' @return The log marginal likelihood (a scalar).
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' lml_brease_mono(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965,
#'                 mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1)
#'
#' @export
lml_brease_mono <- function(y0, y1, N0, N1,
                            mu0 = 0.5, n0 = 2,
                            mue = 0.5, ne = 2) {
  N <- N0 + N1

  lml <- lchoose(N0, y0) + lchoose(N1, y1) -
    lbeta(mu0 * n0, (1 - mu0) * n0) - lbeta(mue * ne, (1 - mue) * ne)

  summands <- lchoose(N1 - y1, 0:(N1 - y1)) +
    lbeta(y0 + y1 + 0:(N1 - y1) + mu0 * n0,
          N - (y0 + y1 + 0:(N1 - y1)) + (1 - mu0) * n0) +
    lbeta(0:(N1 - y1) + mue * ne, y1 + (1 - mue) * ne)

  max_log <- max(summands)
  lml <- lml + log(sum(exp(summands - max_log))) + max_log
  return(lml)
}


#' Log Bayes Factor for the BREASE Model
#'
#' @description Computes the log Bayes factor comparing the sharp null
#'   hypothesis \eqn{H_0: \theta_0 = \theta_1} against the alternative
#'   \eqn{H_1: \theta_0 \neq \theta_1} for the BREASE model.
#'
#' A positive value indicates evidence for H0 (no effect); a negative
#' value indicates evidence for H1 (treatment effect).
#'
#' @inheritParams lml_brease
#'
#' @return The log Bayes factor \eqn{\log BF_{01}} (a scalar). To obtain
#'   the Bayes factor in favor of H1, use \code{exp(-bf_brease(...))}.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial: log BF01
#' lbf01 <- bf_brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
#'                    mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
#'                    mus = 0.3, ns = 1)
#' # BF in favor of treatment effect
#' exp(-lbf01)
#'
#' @export
bf_brease <- function(y0, y1, N0, N1,
                      mu0 = 0.5, n0 = 2,
                      mue = 0.5, ne = 2,
                      mus = 0.5, ns = 2) {
  lml_brease_h0(y0, y1, N0, N1, mu0, n0) -
    lml_brease(y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns)
}


#' Log Bayes Factor for the Monotone BREASE Model
#'
#' @description Computes the log Bayes factor for the monotone (no-harm)
#'   BREASE model.
#'
#' @inheritParams lml_brease
#'
#' @return The log Bayes factor \eqn{\log BF_{01}} (a scalar).
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @export
bf_brease_mono <- function(y0, y1, N0, N1,
                           mu0 = 0.5, n0 = 2,
                           mue = 0.5, ne = 2) {
  lml_brease_h0(y0, y1, N0, N1, mu0, n0) -
    lml_brease_mono(y0, y1, N0, N1, mu0, n0, mue, ne)
}
