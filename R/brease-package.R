#' brease: Causally Sound Priors for Binary Experiments
#'
#' The \code{brease} package implements the BREASE framework for Bayesian
#' analysis of randomized controlled trials with binary treatment and binary
#' outcome, as described in Irons and Cinelli (2025).
#'
#' The BREASE framework parameterizes the likelihood in terms of three
#' clinically meaningful causal quantities: baseline risk (\eqn{\theta_0}),
#' efficacy (\eqn{\eta_e}), and risk of adverse side effects (\eqn{\eta_s}).
#' Independent Beta priors are placed on these parameters, which naturally
#' induces dependence between the treatment and control group risks.
#'
#' The main function is \code{\link{brease}}, which performs posterior
#' inference, computes Bayes factors, and returns an object with
#' print, summary, and plot methods.
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @docType package
#' @name brease-package
#' @importFrom stats rbeta dbeta rbinom density
#' @importFrom stats integrate quantile var cov sd ar
#' @importFrom stats setNames median dnorm rnorm update
#' @importFrom grDevices adjustcolor
#' @importFrom graphics contour plot lines points abline legend text
#'   par axis mtext box image polygon
#' @importFrom graphics filled.contour
"_PACKAGE"
