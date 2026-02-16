#' Exact Posterior Sampling for the BREASE Model
#'
#' @description Draws exact posterior samples from the BREASE model using
#'   data augmentation. This is the default sampling method and requires
#'   no external dependencies.
#'
#' @param n Number of posterior samples to draw.
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
#' @return A list with components:
#' \describe{
#'   \item{theta0}{Posterior samples of baseline risk.}
#'   \item{theta1}{Posterior samples of treatment group risk.}
#'   \item{eta_e}{Posterior samples of efficacy.}
#'   \item{eta_s}{Posterior samples of side-effect risk.}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial data
#' samps <- sample_brease(1000, y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
#'                         mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
#'                         mus = 0.3, ns = 1)
#' hist(samps$theta0, main = "Baseline risk")
#'
#' @export
sample_brease <- function(n, y0, y1, N0, N1,
                          mu0 = 0.5, n0 = 2,
                          mue = 0.5, ne = 2,
                          mus = 0.5, ns = 2) {
  N <- N0 + N1

  # Compute log weights for the data augmentation scheme
  log_weights <- outer(lchoose(y1, 0:y1), lchoose(N1 - y1, 0:(N1 - y1)), "+") +
    lbeta(y1 + y0 + mu0 * n0 + outer(-(0:y1), 0:(N1 - y1), "+"),
          N - y0 - y1 + (1 - mu0) * n0 + outer(0:y1, -(0:(N1 - y1)), "+")) +
    lbeta(rep(1, y1 + 1) %o% 0:(N1 - y1) + mue * ne,
          y1 + (1 - mue) * ne - (0:y1) %o% rep(1, N1 - y1 + 1)) +
    lbeta(0:y1 %o% rep(1, N1 - y1 + 1) + mus * ns,
          N1 - y1 - rep(1, y1 + 1) %o% 0:(N1 - y1) + (1 - mus) * ns)
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)

  # Sample augmented data
  k <- sample(x = 0:(N1 - y1),
              size = n,
              prob = colSums(weights),
              replace = TRUE)
  j <- sapply(seq_len(n), function(i) {
    sample(x = 0:y1, size = 1, prob = weights[, k[i] + 1])
  })

  # Sample parameters from conditional posteriors
  eta_e  <- rbeta(n, k + mue * ne, y1 - j + (1 - mue) * ne)
  eta_s  <- rbeta(n, j + mus * ns, N1 - y1 - k + (1 - mus) * ns)
  theta0 <- rbeta(n, y1 + y0 - j + k + mu0 * n0,
                  N - (y0 + y1 - j + k) + (1 - mu0) * n0)
  theta1 <- theta0 * (1 - eta_e) + (1 - theta0) * eta_s

  list(theta0 = theta0,
       theta1 = theta1,
       eta_e  = eta_e,
       eta_s  = eta_s)
}


#' Exact Posterior Sampling for the Monotone BREASE Model
#'
#' @description Draws exact posterior samples from the monotone (no-harm)
#'   BREASE model, where the side-effect risk \eqn{\eta_s} is constrained
#'   to be zero. This assumes the treatment can only help, not harm.
#'
#' @inheritParams sample_brease
#'
#' @return A list with components:
#' \describe{
#'   \item{theta0}{Posterior samples of baseline risk.}
#'   \item{theta1}{Posterior samples of treatment group risk.}
#'   \item{eta_e}{Posterior samples of efficacy.}
#'   \item{eta_s}{Vector of zeros (no side effects in monotone model).}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' samps <- sample_brease_mono(1000, y0 = 169, y1 = 9, N0 = 20172, N1 = 19965,
#'                              mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1)
#'
#' @export
sample_brease_mono <- function(n, y0, y1, N0, N1,
                               mu0 = 0.5, n0 = 2,
                               mue = 0.5, ne = 2) {
  N <- N0 + N1

  log_weights <- lchoose(N1 - y1, 0:(N1 - y1)) +
    lbeta(y0 + y1 + mu0 * n0 + 0:(N1 - y1),
          N - y0 - y1 + (1 - mu0) * n0 - (0:(N1 - y1))) +
    lbeta(0:(N1 - y1) + mue * ne,
          y1 + (1 - mue) * ne)
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)

  # Sample augmented data
  k <- sample(x = 0:(N1 - y1),
              size = n,
              prob = weights,
              replace = TRUE)

  # Sample parameters
  eta_e  <- rbeta(n, k + mue * ne, y1 + (1 - mue) * ne)
  theta0 <- rbeta(n, y1 + y0 + k + mu0 * n0,
                  N - (y0 + y1 + k) + (1 - mu0) * n0)
  theta1 <- theta0 * (1 - eta_e)

  list(theta0 = theta0,
       theta1 = theta1,
       eta_e  = eta_e,
       eta_s  = rep(0, n))
}


#' Gibbs Sampler for the BREASE Model
#'
#' @description Data-augmented Gibbs sampler for the BREASE model. Useful
#'   as an alternative when the exact sampler is slow (e.g., large sample sizes).
#'
#' @inheritParams sample_brease
#' @param burn_in Number of burn-in iterations to discard (default 1000).
#'
#' @return A data.frame with columns: \code{theta0}, \code{theta1},
#'   \code{eta_e}, \code{eta_s}.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' samps <- gibbs_brease(1000, y0 = 26, y1 = 10, N0 = 11034, N1 = 11037,
#'                        mu0 = 0.5, n0 = 2, mue = 0.3, ne = 1,
#'                        mus = 0.3, ns = 1, burn_in = 500)
#'
#' @export
gibbs_brease <- function(n, y0, y1, N0, N1,
                         mu0 = 0.5, n0 = 2,
                         mue = 0.5, ne = 2,
                         mus = 0.5, ns = 2,
                         burn_in = 1000) {
  N <- N1 + N0
  total_iter <- burn_in + n

  theta0 <- rep(NA, total_iter)
  eta_e  <- rep(NA, total_iter)
  eta_s  <- rep(NA, total_iter)
  theta1 <- rep(NA, total_iter)

  # Initialize
  theta0[1] <- 0.1
  eta_e[1]  <- 0.1
  eta_s[1]  <- 0.1
  theta1[1] <- f_theta1(theta0[1], eta_e[1], eta_s[1])

  for (i in 2:total_iter) {
    # Sample augmented data
    y10 <- rbinom(1, y1, (1 - theta0[i - 1]) * eta_s[i - 1] / theta1[i - 1])
    x11 <- rbinom(1, N1 - y1, theta0[i - 1] * eta_e[i - 1] / (1 - theta1[i - 1]))

    # Sample parameters from full conditionals
    theta0[i] <- rbeta(1, y0 + y1 - y10 + x11 + mu0 * n0,
                       N - (y0 + y1 - y10 + x11) + (1 - mu0) * n0)
    eta_e[i]  <- rbeta(1, x11 + mue * ne, y1 - y10 + (1 - mue) * ne)
    eta_s[i]  <- rbeta(1, y10 + mus * ns, N1 - y1 - x11 + (1 - mus) * ns)
    theta1[i] <- f_theta1(theta0[i], eta_e[i], eta_s[i])
  }

  # Discard burn-in
  idx <- (burn_in + 1):total_iter
  data.frame(theta0 = theta0[idx],
             theta1 = theta1[idx],
             eta_e  = eta_e[idx],
             eta_s  = eta_s[idx])
}
