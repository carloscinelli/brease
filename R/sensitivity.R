#' Bayes Factor Sensitivity Analysis
#'
#' @description Computes the Bayes factor over a grid of prior
#'   hyperparameters to assess sensitivity to the choice of prior.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#' @param mu0 Prior mean for baseline risk (default 0.5).
#' @param n0 Prior sample size for baseline risk (default 2).
#' @param ne Prior sample size for efficacy (default 1).
#' @param ns Prior sample size for side-effect risk (default 1).
#' @param mus_seq Sequence of prior mean values for side-effect risk
#'   (default \code{seq(0.01, 0.99, length = 20)}).
#' @param mue_seq Sequence of prior mean values for efficacy
#'   (default \code{seq(0.01, 0.99, length = 20)}).
#' @param cores Number of cores for parallel computation (default 1,
#'   sequential). Set to a higher value for faster computation.
#'   Requires the \pkg{parallel} package.
#'
#' @return A list with components:
#' \describe{
#'   \item{mus_seq}{The grid of side-effect risk prior means.}
#'   \item{mue_seq}{The grid of efficacy prior means.}
#'   \item{log_bf_seq}{Matrix of log Bayes factors (BF10) at each grid point.}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial sensitivity
#' sens <- bf_seq(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#' contour_bf(sens)
#'
#' @export
bf_seq <- function(y0, y1, N0, N1,
                   mu0 = 0.5, n0 = 2,
                   ne = 1, ns = 1,
                   mus_seq = seq(0.01, 0.99, length = 20),
                   mue_seq = seq(0.01, 0.99, length = 20),
                   cores = 1) {

  f <- function(mue_val, mus_val) {
    -bf_brease(y0, y1, N0, N1, mu0 = mu0, n0 = n0,
               mue = mue_val, ne = ne, mus = mus_val, ns = ns)
  }

  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export necessary functions to workers
    parallel::clusterExport(cl, c("bf_brease", "lml_brease_h0", "lml_brease",
                                   "lml_ib_h0"),
                            envir = environment(bf_brease))

    z <- parallel::parSapply(cl, mue_seq, function(mue_val) {
      sapply(mus_seq, function(mus_val) {
        f(mue_val, mus_val)
      })
    })
  } else {
    z <- sapply(mue_seq, function(mue_val) {
      sapply(mus_seq, function(mus_val) {
        f(mue_val, mus_val)
      })
    })
  }

  out <- list(mus_seq    = mus_seq,
              mue_seq    = mue_seq,
              log_bf_seq = z)
  class(out) <- "bf_seq"
  return(out)
}


#' Contour Plot of Bayes Factor Sensitivity
#'
#' @description Creates a contour plot showing how the Bayes factor varies
#'   across a grid of prior hyperparameters.
#'
#' @param x An object returned by \code{\link{bf_seq}}.
#' @param thr Bayes factor thresholds to highlight (default \code{c(1, 3, 10)}).
#' @param thr_cols Colors for the threshold contours
#'   (default \code{c("red", "orange", "blue")}).
#' @param n_levels Number of contour levels (default 20).
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param cex Text size for contour labels (default 1.5).
#' @param ... Additional arguments passed to \code{contour()}.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' sens <- bf_seq(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#' contour_bf(sens)
#'
#' @export
contour_bf <- function(x,
                       thr = c(1, 3, 10),
                       thr_cols = c("red", "orange", "blue"),
                       n_levels = 20,
                       xlab = "Prior Expected Risk of Side-Effects",
                       ylab = "Prior Expected Efficacy",
                       cex = 1.5, ...) {
  z <- x$log_bf_seq
  mus_seq <- x$mus_seq
  mue_seq <- x$mue_seq

  oldpar <- par(mar = c(4, 4, 1, 1), pty = "s")
  on.exit(par(oldpar))

  levels <- pretty(z, n = n_levels)
  for (j in seq_along(thr)) {
    levels <- levels[!abs(levels - log(thr[j])) < 0.15]
  }

  contour(mus_seq, mue_seq, z, axes = FALSE,
          col = "grey40",
          labels = formatC(exp(levels), 2),
          levels = levels,
          cex = cex,
          xlab = xlab, ylab = ylab, ...)

  for (j in seq_along(thr)) {
    contour(mus_seq, mue_seq, z,
            cex = 2,
            levels = log(thr[j]), label = thr[j],
            col = thr_cols[j], lty = 2, add = TRUE, lwd = 2)
  }

  box()
  axis(1, at = pretty(mus_seq, n = 10), cex.axis = 0.8)
  axis(2, at = pretty(mue_seq, n = 10), cex.axis = 0.8)
}


#' Side-Effect Sensitivity Plot
#'
#' @description Creates a line plot showing how the Bayes factor changes
#'   as a function of the prior expected risk of side effects, for different
#'   prior sample sizes.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#' @param mu0 Prior mean for baseline risk (default 0.5).
#' @param n0 Prior sample size for baseline risk (default 2).
#' @param mue Prior mean for efficacy (default 0.3).
#' @param ne Prior sample size for efficacy (default 1).
#' @param mus_seq Sequence of side-effect risk prior means.
#' @param ns_vals Prior sample size values to compare
#'   (default \code{c(0.5, 1, 2)}).
#' @param lty Line types for each ns value.
#' @param log_scale Logical; if \code{TRUE}, use log scale for y-axis
#'   (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' side_plot(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#'
#' @export
side_plot <- function(y0, y1, N0, N1,
                      mu0 = 0.5, n0 = 2,
                      mue = 0.3, ne = 1,
                      mus_seq = seq(0.001, 0.5, length = 25),
                      ns_vals = c(0.5, 1, 2),
                      lty = c(2, 1, 3),
                      log_scale = FALSE, ...) {

  bf_results <- list()
  for (i in seq_along(ns_vals)) {
    bf_results[[i]] <- sapply(mus_seq, function(mus_val) {
      -bf_brease(y0, y1, N0, N1, mu0 = mu0, n0 = n0,
                 mue = mue, ne = ne, mus = mus_val, ns = ns_vals[i])
    })
  }

  oldpar <- par(mar = c(4, 4, 1, 1), pty = "s")
  on.exit(par(oldpar))

  if (log_scale) {
    f <- function(x) x
    fm1 <- exp
  } else {
    f <- exp
    fm1 <- function(x) x
  }

  range_y <- range(sapply(bf_results, function(x) f(x)))

  plot(mus_seq, f(bf_results[[1]]),
       ylim = range_y * c(1, 1.1),
       type = "n", axes = FALSE,
       xlab = "Prior Expected Risk of Side-Effects",
       ylab = "Evidence against H0",
       cex.lab = 0.8, ...)

  axis(1, at = pretty(mus_seq, n = 10), cex.axis = 0.8)
  yaxis <- pretty(range_y, n = 10)
  axis(2, at = yaxis, labels = formatC(fm1(yaxis), 2), cex.axis = 0.8)

  abline(h = f(log(1)),  col = "red",    lty = 2, lwd = 2)
  abline(h = f(log(3)),  col = "orange", lty = 2, lwd = 2)
  abline(h = f(log(10)), col = "blue",   lty = 2, lwd = 2)

  for (i in seq_along(ns_vals)) {
    lines(mus_seq, f(bf_results[[i]]), type = "l", lty = lty[i])
  }

  legend("topright", xpd = TRUE, bty = "n", ncol = 3,
         title = "Prior sample size", cex = 0.8,
         legend = ns_vals, lty = lty, inset = c(0, 0))
}
