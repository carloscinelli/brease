#' Bayesian Analysis of Binary Experiments Using the BREASE Framework
#'
#' @description Performs Bayesian analysis of a two-arm binary experiment
#'   using the BREASE framework. Parameterizes the model in terms of
#'   baseline risk (\eqn{\theta_0}), efficacy (\eqn{\eta_e}), and adverse
#'   side effects (\eqn{\eta_s}), with treatment group risk given by
#'   \deqn{\theta_1 = \theta_0 (1 - \eta_e) + (1 - \theta_0) \eta_s.}
#'
#'   Independent Beta priors are placed on the causal parameters:
#'   \eqn{\theta_0 \sim \text{Beta}(\mu_0 n_0, (1-\mu_0) n_0)},
#'   \eqn{\eta_e \sim \text{Beta}(\mu_e n_e, (1-\mu_e) n_e)},
#'   \eqn{\eta_s \sim \text{Beta}(\mu_s n_s, (1-\mu_s) n_s)}.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#' @param mu0 Prior mean for baseline risk (default 0.5).
#' @param n0 Prior sample size for baseline risk (default 2).
#' @param mue Prior mean for efficacy (default 0.3).
#' @param ne Prior sample size for efficacy (default 1).
#' @param mus Prior mean for side-effect risk (default 0.3).
#' @param ns Prior sample size for side-effect risk (default 1).
#' @param mono Logical; if \code{TRUE}, use the monotone (no-harm) model
#'   with \eqn{\eta_s = 0} (default \code{FALSE}).
#' @param sampler Character string specifying the sampling method:
#'   \code{"exact"} (default), \code{"gibbs"}, \code{"stan"}, or
#'   \code{"jags"}. The exact and Gibbs samplers are pure R; Stan and
#'   JAGS require the corresponding packages as optional dependencies.
#' @param n_samples Number of posterior samples to draw (default 10000).
#' @param lw Lower quantile for the credible interval (default 0.025).
#' @param up Upper quantile for the credible interval (default 0.975).
#' @param ... Additional arguments passed to the sampler (e.g., \code{burn_in}
#'   for Gibbs/JAGS, \code{iter}, \code{chains} for Stan).
#'
#' @return An object of class \code{"brease"} containing:
#' \describe{
#'   \item{post.samples}{Data.frame of posterior samples for theta0, theta1,
#'     eta_e, eta_s, rr, ef, or, rd.}
#'   \item{post.summaries}{Matrix of posterior mean, median, and credible
#'     interval for each quantity.}
#'   \item{bayes.factors}{List with \code{bf10}, the Bayes factor in favor
#'     of a treatment effect.}
#'   \item{call}{The matched call.}
#'   \item{data}{List with the input data.}
#'   \item{prior}{List with the prior hyperparameters.}
#'   \item{info}{List with method name and settings.}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial (Physicians' Health Study)
#' result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#' result
#' summary(result)
#'
#' # COVID vaccine trial (Pfizer)
#' result <- brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)
#' result
#'
#' # Monotone (no-harm) model
#' result <- brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965, mono = TRUE)
#' result
#'
#' @export
brease <- function(y0, y1, N0, N1,
                   mu0 = 0.5, n0 = 2,
                   mue = 0.3, ne = 1,
                   mus = 0.3, ns = 1,
                   mono = FALSE,
                   sampler = c("exact", "gibbs", "stan", "jags"),
                   n_samples = 10000,
                   lw = 0.025, up = 0.975,
                   ...) {
  # Capture call
  cl <- match.call()
  sampler <- match.arg(sampler)

  # Validate inputs
  check_data(y0, y1, N0, N1)
  check_prior(mu0, n0, mue, ne, mus, ns)

  # Posterior samples
  samps <- switch(sampler,
    exact = {
      if (!mono) {
        sample_brease(n_samples, y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns)
      } else {
        sample_brease_mono(n_samples, y0, y1, N0, N1, mu0, n0, mue, ne)
      }
    },
    gibbs = {
      gibbs_brease(n_samples, y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns, ...)
    },
    stan = {
      stan_brease(n_samples, y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns,
                  mono = mono, ...)
    },
    jags = {
      jags_brease(n_samples, y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns, ...)
    }
  )

  # Convert to data.frame if needed
  if (!is.data.frame(samps)) {
    samps_df <- as.data.frame(samps[c("theta0", "theta1", "eta_e", "eta_s")])
  } else {
    samps_df <- samps
  }

  # Compute treatment effect statistics
  samps_df <- compute_stats(samps_df)

  # Posterior summaries
  post_summaries <- compute_summaries(samps_df, med = 0.5, lw = lw, up = up)

  # Bayes Factors (analytic)
  if (!mono) {
    bf01 <- exp(bf_brease(y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns))
  } else {
    bf01 <- exp(bf_brease_mono(y0, y1, N0, N1, mu0, n0, mue, ne))
  }
  bf10 <- 1 / bf01

  # Build output
  out <- list(
    post.samples   = samps_df,
    post.summaries = post_summaries,
    bayes.factors  = list(bf10 = bf10, bf01 = bf01),
    call           = cl,
    data           = list(y0 = y0, y1 = y1, N0 = N0, N1 = N1),
    prior          = list(mu0 = mu0, n0 = n0, mue = mue, ne = ne,
                          mus = mus, ns = ns),
    info           = list(method = "brease", mono = mono, sampler = sampler)
  )
  class(out) <- "brease"
  return(out)
}


#' Bayesian Analysis Using Independent Beta Priors
#'
#' @description Performs Bayesian analysis of a binary experiment using
#'   independent Beta priors on \eqn{\theta_0} and \eqn{\theta_1}.
#'   This is a standard conjugate approach for comparison with BREASE.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#' @param a Prior shape parameters (vector of length 2 for control and
#'   treatment, default \code{c(1, 1)}).
#' @param b Prior shape parameters (vector of length 2 for control and
#'   treatment, default \code{c(1, 1)}).
#' @param n_samples Number of posterior samples (default 10000).
#' @param lw Lower quantile for the credible interval (default 0.025).
#' @param up Upper quantile for the credible interval (default 0.975).
#'
#' @return An object of class \code{"brease"} (same structure as
#'   \code{\link{brease}}).
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial with uniform (flat) prior
#' result <- brease_ib(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#' result
#'
#' @export
brease_ib <- function(y0, y1, N0, N1,
                      a = c(1, 1), b = c(1, 1),
                      n_samples = 10000,
                      lw = 0.025, up = 0.975) {
  cl <- match.call()
  check_data(y0, y1, N0, N1)

  # Posterior samples (conjugate Beta-Binomial)
  samps <- list(
    theta0 = rbeta(n_samples, y0 + a[1], N0 - y0 + b[1]),
    theta1 = rbeta(n_samples, y1 + a[2], N1 - y1 + b[2])
  )
  samps_df <- as.data.frame(samps)

  # Compute treatment effect statistics
  samps_df <- compute_stats(samps_df)

  # Posterior summaries
  post_summaries <- compute_summaries(samps_df, med = 0.5, lw = lw, up = up)

  # Bayes Factor (analytic)
  bf10 <- exp(bf_ib_analytic(y0, y1, N0, N1, a = a[1], b = b[1]))

  out <- list(
    post.samples   = samps_df,
    post.summaries = post_summaries,
    bayes.factors  = list(bf10 = bf10, bf01 = 1 / bf10),
    call           = cl,
    data           = list(y0 = y0, y1 = y1, N0 = N0, N1 = N1),
    prior          = list(a = a, b = b),
    info           = list(method = "ib", mono = FALSE, sampler = "conjugate")
  )
  class(out) <- "brease"
  return(out)
}


#' Bayesian Analysis Using the Logit-Transform Prior
#'
#' @description Performs Bayesian analysis of a binary experiment using
#'   the logit-transform (LT) prior, which places normal priors on the
#'   log-odds baseline and log-odds ratio. Requires the \pkg{abtest}
#'   package.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#' @param mu_psi Prior mean for the log-odds ratio (default 0).
#' @param sigma_psi Prior standard deviation for the log-odds ratio (default 1).
#' @param mu_beta Prior mean for the log-odds baseline (default 0).
#' @param sigma_beta Prior standard deviation for the log-odds baseline (default 1).
#' @param n_samples Number of posterior samples (default 10000).
#' @param lw Lower quantile for the credible interval (default 0.025).
#' @param up Upper quantile for the credible interval (default 0.975).
#'
#' @return An object of class \code{"brease"} (same structure as
#'   \code{\link{brease}}).
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' \dontrun{
#' # Requires the abtest package
#' result <- brease_lt(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#' result
#' }
#'
#' @export
brease_lt <- function(y0, y1, N0, N1,
                      mu_psi = 0, sigma_psi = 1,
                      mu_beta = 0, sigma_beta = 1,
                      n_samples = 10000,
                      lw = 0.025, up = 0.975) {
  cl <- match.call()
  check_data(y0, y1, N0, N1)

  if (!requireNamespace("abtest", quietly = TRUE)) {
    stop("Package 'abtest' is required for brease_lt(). ",
         "Please install it with: install.packages('abtest')",
         call. = FALSE)
  }

  # Use abtest to compute posterior
  ab <- abtest::ab_test(
    data = data.frame(y1 = y0, n1 = N0, y2 = y1, n2 = N1),
    nsamples = n_samples,
    prior_par = list(mu_psi = mu_psi, sigma_psi = sigma_psi,
                     mu_beta = mu_beta, sigma_beta = sigma_beta)
  )
  samps <- abtest::get_post_samples(ab)
  samps <- samps[c("p1", "p2")]
  colnames(samps) <- c("theta0", "theta1")

  # Compute treatment effect statistics
  samps_df <- compute_stats(as.data.frame(samps))

  # Posterior summaries
  post_summaries <- compute_summaries(samps_df, med = 0.5, lw = lw, up = up)

  # Bayes Factor
  bf10 <- exp(ab$logbf$bf10)

  out <- list(
    post.samples   = samps_df,
    post.summaries = post_summaries,
    bayes.factors  = list(bf10 = bf10, bf01 = 1 / bf10),
    call           = cl,
    data           = list(y0 = y0, y1 = y1, N0 = N0, N1 = N1),
    prior          = list(mu_psi = mu_psi, sigma_psi = sigma_psi,
                          mu_beta = mu_beta, sigma_beta = sigma_beta),
    info           = list(method = "lt", mono = FALSE, sampler = "abtest")
  )
  class(out) <- "brease"
  return(out)
}
