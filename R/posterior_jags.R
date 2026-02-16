#' JAGS MCMC Sampling for the BREASE Model
#'
#' @description Uses JAGS (via the \pkg{rjags} package) to draw posterior
#'   samples from the BREASE model. Requires \pkg{rjags} to be installed.
#'
#' @inheritParams sample_brease
#' @param burn_in Number of burn-in iterations (default 1000).
#' @param ... Additional arguments passed to \code{rjags::jags.model()}.
#'
#' @return A data.frame with columns: \code{theta0}, \code{theta1},
#'   \code{eta_e}, \code{eta_s}.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @export
jags_brease <- function(n = 10000, y0, y1, N0, N1,
                        mu0 = 0.5, n0 = 2,
                        mue = 0.5, ne = 2,
                        mus = 0.5, ns = 2,
                        burn_in = 1000, ...) {

  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' is required for sampler = 'jags'. ",
         "Please install it with: install.packages('rjags')",
         call. = FALSE)
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for sampler = 'jags'. ",
         "Please install it with: install.packages('coda')",
         call. = FALSE)
  }

  a0 <- mu0 * n0
  b0 <- (1 - mu0) * n0
  ae <- ne * mue
  be <- ne * (1 - mue)
  as_ <- ns * mus
  bs <- ns * (1 - mus)

  model_code <- "model{
    y1 ~ dbinom(theta1, N1);
    y0 ~ dbinom(theta0, N0);
    theta0  ~ dbeta(a0, b0);
    etae    ~ dbeta(ae, be);
    etas    ~ dbeta(as_, bs);
    theta1 <- theta0 * (1 - etae) + etas * (1 - theta0);
  }"

  jags_data <- list(y1 = y1, y0 = y0, N1 = N1, N0 = N0,
                    a0 = a0, b0 = b0, ae = ae, be = be,
                    as_ = as_, bs = bs)

  model <- rjags::jags.model(
    file = textConnection(model_code),
    data = jags_data,
    quiet = TRUE, ...
  )
  stats::update(model, burn_in)
  samps <- rjags::coda.samples(
    model = model,
    variable.names = c("theta1", "theta0", "etae", "etas"),
    n.iter = n
  )
  samps <- data.frame(samps[[1]])

  data.frame(theta0 = samps$theta0,
             theta1 = samps$theta1,
             eta_e  = samps$etae,
             eta_s  = samps$etas)
}
