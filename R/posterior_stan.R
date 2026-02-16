#' Stan MCMC Sampling for the BREASE Model
#'
#' @description Uses Stan (via the \pkg{rstan} package) to draw posterior
#'   samples from the BREASE model. Requires \pkg{rstan} to be installed.
#'
#' @inheritParams sample_brease
#' @param iter Total number of iterations per chain (default 10000).
#' @param chains Number of MCMC chains (default 4).
#' @param warmup Number of warmup iterations per chain (default \code{iter / 2}).
#' @param mono Logical; if \code{TRUE}, use the monotone (no-harm) model
#'   (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{rstan::sampling()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{theta0}{Posterior samples of baseline risk.}
#'   \item{theta1}{Posterior samples of treatment group risk.}
#'   \item{eta_e}{Posterior samples of efficacy.}
#'   \item{eta_s}{Posterior samples of side-effect risk (0 if monotone).}
#'   \item{stanfit}{The raw stanfit object for diagnostics.}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @export
stan_brease <- function(n = 10000, y0, y1, N0, N1,
                        mu0 = 0.5, n0 = 2,
                        mue = 0.5, ne = 2,
                        mus = 0.5, ns = 2,
                        iter = 10000, chains = 4, warmup = iter %/% 2,
                        mono = FALSE, ...) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for sampler = 'stan'. ",
         "Please install it with: install.packages('rstan')",
         call. = FALSE)
  }

  stan_data <- list(N1 = as.integer(N1), y1 = as.integer(y1),
                    N0 = as.integer(N0), y0 = as.integer(y0),
                    mu0 = mu0, n0 = n0,
                    mue = mue, ne = ne,
                    mus = mus, ns = ns)

  if (mono) {
    code <- .mono_stan_code
  } else {
    code <- .brease_stan_code
  }

  model <- rstan::stan_model(model_code = code)
  fit <- rstan::sampling(model, data = stan_data,
                         iter = iter, chains = chains, warmup = warmup, ...)
  samples <- rstan::extract(fit)

  out <- list(theta0 = as.numeric(samples$theta0),
              theta1 = as.numeric(samples$theta1),
              eta_e  = as.numeric(samples$etae))
  if (!mono) {
    out$eta_s <- as.numeric(samples$etas)
  } else {
    out$eta_s <- rep(0, length(out$theta0))
  }
  out$stanfit <- fit
  out
}


# Stan model code (full BREASE model)
.brease_stan_code <- '
data {
  int<lower=1> N1;
  int<lower=1> N0;
  int<lower=0,upper=N1> y1;
  int<lower=0,upper=N0> y0;
  real<lower=0,upper=1> mu0;
  real<lower=0,upper=1> mue;
  real<lower=0,upper=1> mus;
  real<lower=0> n0;
  real<lower=0> ne;
  real<lower=0> ns;
}
parameters {
  real<lower=0,upper=1> theta0;
  real<lower=0,upper=1> etae;
  real<lower=0,upper=1> etas;
}
transformed parameters {
  real<lower=0,upper=1> theta1 = (1-etae)*theta0 + etas*(1-theta0);
}
model {
  target += beta_lpdf(theta0 | mu0*n0, (1-mu0)*n0);
  target += beta_lpdf(etae | mue*ne, (1-mue)*ne);
  target += beta_lpdf(etas | mus*ns, (1-mus)*ns);
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'

# Stan model code (monotone/no-harm model)
.mono_stan_code <- '
data {
  int<lower=1> N1;
  int<lower=1> N0;
  int<lower=0,upper=N1> y1;
  int<lower=0,upper=N0> y0;
  real<lower=0,upper=1> mu0;
  real<lower=0,upper=1> mue;
  real<lower=0,upper=1> mus;
  real<lower=0> n0;
  real<lower=0> ne;
  real<lower=0> ns;
}
parameters {
  real<lower=0,upper=1> theta0;
  real<lower=0,upper=1> etae;
}
transformed parameters {
  real<lower=0,upper=1> theta1 = (1-etae)*theta0;
}
model {
  target += beta_lpdf(theta0 | mu0*n0, (1-mu0)*n0);
  target += beta_lpdf(etae | mue*ne, (1-mue)*ne);
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'
