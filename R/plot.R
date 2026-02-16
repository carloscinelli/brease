#' Plot Method for BREASE Objects
#'
#' @description Creates density plots of posterior distributions from a
#'   BREASE analysis.
#'
#' @param x An object of class \code{"brease"}.
#' @param type Character string specifying the type of plot:
#'   \code{"posterior"} (default) for posterior densities, or
#'   \code{"prior"} for prior predictive densities.
#' @param which Character vector specifying which quantities to plot
#'   (default \code{c("theta0", "theta1")}). Can include any column
#'   from the posterior samples.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
plot.brease <- function(x, type = c("posterior", "prior"),
                        which = c("theta0", "theta1"), ...) {
  type <- match.arg(type)

  if (type == "posterior") {
    samps <- x$post.samples
  } else {
    # Simulate from the prior
    method <- x$info$method
    if (method == "brease") {
      p <- x$prior
      samps <- sim_brease(10000, mu0 = p$mu0, n0 = p$n0,
                           mue = p$mue, ne = p$ne,
                           mus = p$mus, ns = p$ns)
      samps <- compute_stats(samps)
    } else if (method == "ib") {
      p <- x$prior
      samps <- sim_ib(10000, mu0 = p$a[1] / (p$a[1] + p$b[1]),
                        n0 = p$a[1] + p$b[1],
                        mu1 = p$a[2] / (p$a[2] + p$b[2]),
                        n1 = p$a[2] + p$b[2])
      samps <- compute_stats(samps)
    } else if (method == "lt") {
      p <- x$prior
      samps <- sim_lt(10000, mu_psi = p$mu_psi, var_psi = p$sigma_psi^2,
                        mu_beta = p$mu_beta, var_beta = p$sigma_beta^2)
      samps <- compute_stats(samps)
    }
  }

  # Filter to requested quantities
  which <- which[which %in% names(samps)]
  n_plots <- length(which)
  if (n_plots == 0) {
    warning("No valid quantities to plot.")
    return(invisible(x))
  }

  # Labels
  labels <- c(theta0 = expression(theta[0]),
              theta1 = expression(theta[1]),
              eta_e  = expression(eta[e]),
              eta_s  = expression(eta[s]),
              rr     = "Risk Ratio",
              ef     = "Efficacy Fraction",
              or     = "Odds Ratio",
              rd     = "Risk Difference")

  oldpar <- par(mfrow = c(1, n_plots),
                mar = c(4, 4, 2.5, 1))
  on.exit(par(oldpar))

  title_prefix <- if (type == "posterior") "Posterior" else "Prior"

  for (w in which) {
    vals <- samps[[w]]
    d <- density(vals, na.rm = TRUE)
    lab <- if (w %in% names(labels)) labels[[w]] else w
    plot(d, main = paste(title_prefix, "-", w),
         xlab = lab, ylab = "Density",
         col = "steelblue", lwd = 2,
         font.main = 1, ...)
    polygon(d, col = adjustcolor("steelblue", alpha.f = 0.3), border = NA)
  }

  invisible(x)
}
